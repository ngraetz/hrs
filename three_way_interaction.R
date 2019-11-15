library(data.table)
library(ggplot2)
library(formatR)
library(dplyr)
library(gridExtra)
library(lme4)
library(survey)
library(mice)
library(survey)
library(purrr)
library(stargazer)
library(sjPlot)

repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
hrs <- readRDS('hrs_imputed_v2.RDS')

## Restrict race
hrs <- hrs[race %in% c('white','black'), ]
hrs[, race := factor(race, levels=c('white','black'))]
hrs[, edu_cat := factor(edu_cat, levels=c('highschool','less_highschool','some_college','college'))]

outcome_var <- 'cognitive'
REML <- FALSE
weight <- FALSE

numeric_vars <- c('baseline_wealth','baseline_income') 
factor_vars <- c('as.factor(female)','as.factor(cohort_group)','edu_cat')

## MATCHING
confounders <- c('age', 'female', 'cohort_group_1930', 'cohort_group_1940', 'baseline_cog', 'edu_cat_less_highschool','edu_cat_highschool','edu_cat_some_college', 'baseline_wealth', 'baseline_income')
treatment <- 'race_black'

## Calculate baseline cognition and only keep those where we observe cognition at "baseline" (before age 60)
model_hrs <- copy(hrs)
model_hrs[, age_obs := age]
model_hrs[, age_obs := min(age_obs, na.rm = T), by = id]
model_hrs[age==age_obs, baseline_cog := get(outcome_var)]
model_hrs[, baseline_cog := max(baseline_cog, na.rm = T), by = id] ## Repeat within individual
model_hrs <- model_hrs[age_obs < 60, ]
model_hrs[baseline_cog<15, baseline := 'low']
model_hrs[baseline_cog %in% 15:19, baseline := 'middle']
model_hrs[baseline_cog %in% 20:35, baseline := 'high']
model_hrs <- model_hrs[!is.infinite(baseline_cog), ]
model_hrs[, original_baseline_cog := baseline_cog]

## Try calculating baseline wealth and log income as well
model_hrs[age==age_obs, baseline_wealth := wealth]
model_hrs[, baseline_wealth := max(baseline_wealth, na.rm = T), by = id] ## Repeat within individual
model_hrs <- model_hrs[!is.infinite(baseline_wealth), ]
model_hrs[, original_baseline_wealth := baseline_wealth]

model_hrs[age==age_obs, baseline_income := log_income]
model_hrs[, baseline_income := max(baseline_income, na.rm = T), by = id] ## Repeat within individual
model_hrs <- model_hrs[!is.infinite(baseline_income), ]
model_hrs[, original_baseline_income := baseline_income]

## Subset to respondents with at least 3 responses on outcome.
id_nums <- model_hrs[, .N , by = id]
dim(id_nums[N == 1, ])
model_hrs <- model_hrs[id %in% id_nums[N > 3, id], ]
# model_hrs[, id_factor := as.factor(id)]

## Rescale all variables to mean/sd=0/1. Keep track of means/sds to transform back after predicting.
scale_vars <- c(outcome_var, numeric_vars, 'baseline_cog')
scale_vars <- scale_vars[!is.na(scale_vars)]
means <- model_hrs[, lapply(.SD, mean, na.rm=TRUE), .SDcols=scale_vars]
means <- melt(means, measure.vars=scale_vars)
sds <- model_hrs[, lapply(.SD, sd, na.rm=TRUE), .SDcols=scale_vars]
sds <- melt(sds, measure.vars=scale_vars)
setnames(sds, 'value', 'sd')
setnames(means, 'value', 'mean')
scales <- merge(means,sds)
model_hrs <- as.data.table(model_hrs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=scale_vars)) # Subtract mean, divide by SD.

## Center linear age at 50 to interpret baseline intercept
model_hrs[, center_age := (age - 50) / 5]

## If we want age polynomials, we need to calculate orthogonal polynomials. Because our ages are relatively high, age and age^2
## are basically perfectly correlated so the model can't converge.
age_polys <- poly(model_hrs[, center_age], degree = 2)
model_hrs[, age_poly_1 := age_polys[,1]]
model_hrs[, age_poly_2 := age_polys[,2]]
model_hrs <- as.data.table(model_hrs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=c('age_poly_1','age_poly_2'))) # Subtract mean, divide by SD.

## Fit model, pull out random effects to calculate individual trajectories. Coef() gives random + fixed effects.
model_hrs[, race := factor(race, levels=c('white','black','hispanic','other'))]
model_hrs[, id_factor := as.factor(id)]
model_hrs <- model_hrs[!is.na(pweight) & pweight > 0, ]
model_hrs[, pweight := pweight / 10000]

formula1 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars,numeric_vars), collapse=' + '), ' + center_age*baseline_cog*race + (center_age||id_factor)'))
model1 <- lmer(formula1, weights=pweight, data=model_hrs, REML=FALSE)

formula2 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars,numeric_vars), collapse=' + '), ' + baseline_cog + center_age*edu_cat*race + (center_age||id_factor)'))
model2 <- lmer(formula2, weights=pweight, data=model_hrs, REML=FALSE)

## TRY MATCHING BASELINE OBSERVATIONS 
source('C:/Users/ngraetz/Documents/repos/hrs/causal_functions.R')
d <- model_hrs[age==age_obs & race %in% c('black','white') & !is.na(edu_years), ]
for(var in c('cohort_group','race','edu_cat')) {
  for(v in unique(d[, get(var)])) {
    d[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}

## Use propensity score matching with Mahalanobis distance and caliper.
ps_formula <- as.formula(paste0(treatment, ' ~ ', paste(confounders, collapse = ' + ')))
ps_mod <- glm(ps_formula,
              family=binomial,
              x=TRUE,
              y=TRUE,
              data=d)
d[, prop_score := predict(ps_mod)]
max.ps.control <- max(d[get(treatment)==F, prop_score])
min.ps.treated <- min(d[get(treatment)==T, prop_score])
d[get(treatment)==F & prop_score < min.ps.treated, drop := 1]
d[get(treatment)==T & prop_score > max.ps.control, drop := 1]
message(paste0('Dropping ', length(d[drop==1, drop]), ' obs with outlier propensity scores.'))
d <- d[is.na(drop), ]

## Set up for matching.
# Matrix of confounders.
Xmat <- d[, confounders, with=F]
# Matrix of covariates to include in the Mahalanobis distance (just use all).
Xmatmahal <- d[, confounders, with=F]
# Rank based Mahalanobis distance
distmat <- smahal(d[, get(treatment)], Xmatmahal)
# Add caliper
distmat2 <- addcaliper(distmat, d[, get(treatment)], d[, prop_score], calipersd=.5)
# Name the rows and columns of distance matrix by the subject numbers in treated
d[, id := 1:dim(d)[1]]
rownames(distmat2)=d[get(treatment)==T, id]
colnames(distmat2)=d[get(treatment)==F, id]

## Calculate matches (1, 2, 3, full)
options("optmatch_max_problem_size" = Inf)
# lapply(c(1,2,3), match_controls)
match_controls(1)

## Check for acceptable balance
make_diffs <- function(n_controls) {
  diffs.after <- rbindlist(lapply(confounders, calc_stand_diff,
                                  treatment = treatment,
                                  data = d[!is.na(get(paste0('matchvec_', n_controls))), ]))
  diffs.after[, matching := paste0(n_controls, ' matched controls')]
  return(diffs.after)
}
all.diffs <- make_diffs('1')
diffs.before <- rbindlist(lapply(confounders, calc_stand_diff, treatment = treatment, data = d))
diffs.before[, matching := 'Before matching']
all.diffs <- rbind(all.diffs, diffs.before)

## Make Love Plot
gg.love <- ggplot() + 
  geom_jitter(data=all.diffs,
              aes(x=diff,
                  y=variable,
                  fill=matching),
              width = 0.1,
              height = 0.1,
              size=5,
              shape=21) + 
  geom_vline(xintercept = .2, linetype='dashed') + 
  geom_vline(xintercept = -.2, linetype='dashed') + 
  geom_vline(xintercept = 0) + 
  labs(y='',x='Standardized difference',
       title='Love plot of standardized differences for all confounders\nMatching via rank-based Mahalanobis distance with a propensity score caliper', fill='') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

## Merge matches back on to full dataset and examine raw trajectories over age.
d2 <- merge(model_hrs, d[!is.na(matchvec_1), c('id_factor','matchvec_1')], by='id_factor')
d2[, unique_ids_race := uniqueN(id_factor), by = 'race']
d2[, unique_ids_race_edu := uniqueN(id_factor), by = c('race','edu_cat')]
trends <- d2[, list(cognitive=mean(cognitive, na.rm=T)), by=c('age','race')]
trends <- merge(trends, unique(d2[,c('race','unique_ids_race')]), by='race')
gg.race <- ggplot(data=trends,
                  aes(x=age,
                      y=cognitive)) + 
            geom_line(aes(color=race,
                          fill=race),
                      size=2) + 
            geom_smooth(aes(color=race,
                            fill=race),
                        alpha=0.5) + 
            geom_label(data = trends[age==60, ],
                       aes(x = age,
                           y = cognitive,
                           label = unique_ids_race,
                           group = race)) + 
            theme_bw()
trends <- d2[, list(cognitive=mean(cognitive, na.rm=T)), by=c('age','race','edu_cat')]
trends <- trends[edu_cat %in% c('college','less_highschool')]
trends[, race_edu := paste0(race,'_',edu_cat)]
trends <- merge(trends, unique(d2[,c('race','edu_cat','unique_ids_race_edu')]), by=c('race','edu_cat'))
gg.race.edu <- ggplot(data=trends,
                  aes(x=age,
                      y=cognitive)) + 
  geom_line(aes(color=race_edu,
                fill=race_edu),
            size=2) + 
  geom_smooth(aes(color=race_edu,
                  fill=race_edu),
              alpha=0.5) + 
  geom_label(data = trends[age==60, ],
             aes(x = age,
                 y = cognitive,
                 label = unique_ids_race_edu,
                 group = race_edu)) + 
  scale_color_manual(values = c('#3182bd','#9ecae1','#de2d26','#fc9272')) + 
  scale_fill_manual(values = c('#3182bd','#9ecae1','#de2d26','#fc9272')) + 
  theme_bw()

## Make all plots
pdf('matched_pairs_edu_cat.pdf', width=12, height=8)
print(gg.love)
print(gg.race)
print(gg.race.edu)
dev.off()

## Try modelling with fixed effect on matched pairs.
formula3 <- as.formula(paste0(outcome_var, ' ~ matchvec_1 + center_age*edu_cat*race + (center_age||id_factor)'))
model3 <- lmer(formula3, weights=pweight, data=d2, REML=FALSE)

## Make model table
tab_model(list(model1,model2,model3), dv.labels = c('Model 1','Model 2','Matched'), show.ci=FALSE, p.style='asterisk', file='three_way_interaction_edu_cat.html')

## Try making predicted scores
coef_table <- as.data.table(expand.grid(sort(unique(model_hrs[, center_age])), c('white','black'), c('less_highschool','highschool','some_college','college')))
setnames(coef_table, c('center_age','race','edu_cat'))
for(e in unique(coef_table[, edu_cat])) coef_table[, (paste0('edu_cat',e)) := ifelse(edu_cat==e, 1, 0)]
coef_table[, raceblack := ifelse(race=='black',1,0)]
for(e in c('less_highschool','highschool','some_college','college')) {
  ## EDU * AGE DUMMIES
  coef_table[, (paste0('center_age:edu_cat',e)) := center_age * get(paste0('edu_cat',e))]
  ## EDU * RACE DUMMIES
  coef_table[, (paste0('edu_cat',e,':raceblack')) := get(paste0('edu_cat',e)) * raceblack]
  ## RACE * AGE DUMMIES
  coef_table[, ('center_age:raceblack') := center_age * raceblack]
  ## EDU * AGE * RACE DUMMIES
  coef_table[, (paste0('center_age:edu_cat',e, ':raceblack')) := center_age * get(paste0('edu_cat',e)) * raceblack]
}
coef_table[, ("(Intercept)") := 1]
betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway.
betas <- betas[, !c('matchvec_1'), with=F]
beta_names <- names(betas)
for(d in 1:1000) {
  betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
  names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway.
  betas <- betas[, !c('matchvec_1'), with=F]
  beta_names <- names(betas)
  template <- coef_table[, beta_names, with=F]
  setcolorder(template, beta_names)
  coef_table[, (paste0('draw',d)) := as.numeric(as.matrix(betas) %*% t(as.matrix(template)))]
}
coef_table[, pred_mean := apply(.SD,1,mean), .SDcols=paste0('draw',1:1000)]
coef_table[, pred_lower := apply(.SD,1,quantile,probs=0.025), .SDcols=paste0('draw',1:1000)]
coef_table[, pred_upper := apply(.SD,1,quantile,probs=0.975), .SDcols=paste0('draw',1:1000)]
coef_table[, race_edu := paste0(race,'_',edu_cat)]
coef_table[, center_age := center_age * 5 + 50]
coef_table[race=='black', race := 'Black']
coef_table[race=='white', race := 'White']
coef_table[edu_cat=='college', edu_cat := 'College']
coef_table[edu_cat=='some_college', edu_cat := 'Some College']
coef_table[edu_cat=='highschool', edu_cat := 'High School']
coef_table[edu_cat=='less_highschool', edu_cat := 'Less High School']
pdf('C:/Users/ngraetz/Documents/repos/hrs/predicted_means.pdf', width=10, heigh=8)
ggplot(data=coef_table, aes(x=center_age,y=pred_mean)) + 
  geom_ribbon(aes(ymin=pred_lower,
                  ymax=pred_upper,
                  fill=race),
              color='black',
              alpha=0.6,
              size=0.2) + 
  geom_line(aes(group=race),color='black', size=1) + 
  facet_wrap(~edu_cat) + 
  theme_bw() + 
  labs(y='Predicted mean cognition score',x='Age',fill='Race') + 
  theme(strip.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.text = element_text(size = 20),
        legend.key.size = unit(3,'line'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
dev.off()

## Make coefficient comparison graph
betas <- as.data.table(MASS::mvrnorm(1000, mu = fixef(model3), Sigma = vcov(model3)))
for(e in paste0('edu_cat',c('college','some_college','highschool','less_highschool'))) {
  for(r in c('white','black')) {
    ## AGE*EDU
    if(e!='edu_cathighschool') betas[, age_edu := get(paste0('center_age:',e))]
    if(e=='edu_cathighschool') betas[, age_edu := 0]
    ## AGE*RACE
    if(r=='black') betas[, age_race := get('center_age:raceblack')]
    if(r!='black') betas[, age_race := 0]
    ## AGE*RACE*EDU
    if(r=='black' & e!='edu_cathighschool') betas[, age_race_edu := get(paste0('center_age:',e,':raceblack'))]
    if(r!='black') betas[, age_race_edu := 0]
    ## FULL COEFFICIENT
    betas[, (paste0(r,'_',e)) := center_age + age_edu + age_race + age_race_edu]
  }
}
betas <- betas[, c(paste0('white_',paste0('edu_cat',c('college','some_college','highschool','less_highschool'))),
                   paste0('black_',paste0('edu_cat',c('college','some_college','highschool','less_highschool'))))]
betas <- melt(betas, value.name = 'coef')
betas[, draw := seq(1:.N), by='variable']
betas[, mean := mean(coef), by='variable']
betas[, lower := quantile(coef,probs=0.025), by='variable']
betas[, upper := quantile(coef,probs=0.975), by='variable']
betas <- unique(betas[, c('variable','mean','upper','lower')])
betas[, variable := gsub('edu_cat','',variable)]
betas[, race := tstrsplit(variable,'_',keep=1)]
betas[, edu := tstrsplit(variable,'_',keep=2)]
betas[edu=='some', edu := 'Some College']
betas[edu=='less', edu := 'Less High School']
betas[edu=='college', edu := 'College']
betas[edu=='highschool', edu := 'High School']
betas[, edu := factor(edu, levels=c('Less High School','High School','Some College','College'))]
pdf('C:/Users/ngraetz/Documents/repos/hrs/predicted_coefficients.pdf', width=8, heigh=6)
ggplot(data=betas, aes(x=edu,y=mean)) +
  geom_errorbar(aes(ymin=lower,ymax=upper,group=race), size=1) + 
  geom_point(aes(fill=race), size=8, shape=21) +
  labs(y='Fitted coefficient on age (5-year increase)',x='',fill='Race') + 
  theme_bw()
dev.off()

## Descriptive table
table1 <- copy(model_hrs)
for(var in c('cohort_group','race','edu_cat')) {
  for(v in unique(table1[, get(var)])) {
    table1[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}
measures <- c('female','race_white','race_black','race_hispanic','race_other','cohort_group_1930','cohort_group_1940','cohort_group_1950','original_baseline_cog','edu_cat_less_highschool','edu_cat_highschool','edu_cat_some_college','edu_cat_college')
table1[, n := 1]
totals <- table1[!is.na(age_group), list(N=sum(n)), by='age_group']
table1 <- table1[!is.na(age_group), lapply(.SD, mean, na.rm=TRUE), by='age_group', .SDcols=measures]
table1 <- merge(table1, totals, by='age_group')
table1 <- melt(table1, id.vars = 'age_group')
table1[!(variable %in% c('original_baseline_cog','N')), value := value * 100]
table1[, value := ifelse(variable=='N', round(value), round(value,1))]
table1[, value := as.character(value)]
table1[!grepl('[.]',value) & value!='0' & variable!='N', value := paste0(value, '.0')]
table1 <- dcast(table1, variable ~ age_group, value.var='value')
write.csv(table1, 'table1.csv', row.names = F)
saveRDS(table1, 'table1.rds')

table1 <- copy(model_hrs)
for(var in c('cohort_group','race','edu_cat')) {
  for(v in unique(table1[, get(var)])) {
    table1[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}
table1 <- table1[!is.na(age_group), lapply(.SD, weighted.mean, w=pweight, na.rm=TRUE), by='age_group', .SDcols=measures]
table1 <- merge(table1, totals, by='age_group')
table1 <- melt(table1, id.vars = 'age_group')
table1[!(variable %in% c('original_baseline_cog','N')), value := value * 100]
table1[, value := ifelse(variable=='N', round(value), round(value,1))]
table1[, value := as.character(value)]
table1[!grepl('[.]',value) & value!='0' & variable!='N', value := paste0(value, '.0')]
table1 <- dcast(table1, variable ~ age_group, value.var='value')
write.csv(table1, 'table1_weighted.csv', row.names = F)
saveRDS(table1, 'table1_weighted.rds')
