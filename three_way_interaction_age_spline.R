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
library(lmerTest)
library(MuMIn)
library(flextable)

repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
hrs <- readRDS('hrs_imputed_v3.RDS')
oldest <- FALSE
drop_oldest <- FALSE
file_tag <- 'age_interaction'

## Restrict race
hrs <- hrs[race %in% c('white','black'), ]
hrs[, race := factor(race, levels=c('white','black'))]
hrs[, edu_cat := factor(edu_cat, levels=c('highschool','less_highschool','some_college','college'))]

outcome_var <- 'cognitive'
REML <- FALSE
weight <- FALSE

# 'baseline_wealth','baseline_income','baseline_cognition'
numeric_vars <- c('edu_years','baseline_wealth','baseline_income','baseline_cog') 
factor_vars <- c('as.factor(female)','as.factor(cohort_group)')

## MATCHING
ifelse(oldest, confounders <- c('age', 'female', 'cohort_group_1930', 'baseline_cog', 'edu_years', 'baseline_wealth', 'baseline_income'), confounders <- c('age', 'female', 'cohort_group_1930', 'cohort_group_1940', 'baseline_cog', 'edu_years', 'baseline_wealth', 'baseline_income'))
treatment <- 'race_black'

## Calculate baseline cognition and only keep those where we observe cognition at "baseline" (before age 60)
model_hrs <- copy(hrs)
if(oldest) model_hrs <- model_hrs[age>=70]
if(drop_oldest) model_hrs <- model_hrs[age<70]
model_hrs[, age_obs := age]
model_hrs[, age_obs := min(age_obs, na.rm = T), by = id]
model_hrs[age==age_obs, baseline_cog := get(outcome_var)]
model_hrs[, baseline_cog := max(baseline_cog, na.rm = T), by = id] ## Repeat within individual
ifelse(oldest, model_hrs <- model_hrs[age_obs < 75, ], model_hrs <- model_hrs[age_obs < 60, ])
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
model_hrs[, original_cognitive := cognitive]
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
ifelse(oldest, model_hrs[, center_age := (age - 70) / 5], model_hrs[, center_age := (age - 50) / 5])

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

## Term for practice effects
model_hrs[, practice := ifelse(age==age_obs, 0, 1)]

## Term for over/under age 70
model_hrs[, over_70 := ifelse(center_age>=4, 1, 0)]

# formula1 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars,numeric_vars), collapse=' + '), ' + center_age*baseline_cog*race + (center_age||id_factor)'))
# model1 <- lmerTest::lmer(formula1, weights=pweight, data=model_hrs, REML=FALSE)

# M1: basic controls (age, sex, race, cohort)
# M2: basic controls + baseline cognition
# M3: M2 + early-life SES, early-life SES*age
# M4: M2 + adult SES + race*edu, race*edu*age
# M5: fully adjusted
# formula2 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars,numeric_vars), collapse=' + '), ' + center_age*edu_years*race + (center_age||id_factor)'))

for(v in c('cognitive','baseline_cog','baseline_wealth','baseline_income','edu_years_mother','edu_years','cohort_group')) model_hrs <- model_hrs[!is.na(get(v)), ]
model_hrs[, demean_cognitive := mean(cognitive), by='id']
model_hrs[, demean_cognitive := cognitive - demean_cognitive]

formula1 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
  bs(center_age, knots=4, degree=1)*race + 
  (center_age||id_factor)
formula2 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + bs(center_age, knots=4, degree=1)*race + 
  baseline_cog + (center_age||id_factor)
# formula3 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   center_age*edu_years_mother*race + 
#   baseline_cog + (center_age||id_factor)
# formula3 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   center_age*over_70*race +
#   center_age*edu_years + 
#   baseline_cog + baseline_wealth + baseline_income + (center_age||id_factor)
formula3 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
  bs(center_age, knots=4, degree=1)*race +
  bs(center_age, knots=4, degree=1)*edu_years + 
  baseline_cog + baseline_wealth + baseline_income + (center_age||id_factor)
# formula5 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   center_age*edu_years*race + center_age*edu_years_mother*race + 
#   baseline_cog + baseline_wealth + baseline_income + (center_age||id_factor)

# ols1 <- lm(demean_cognitive ~ as.factor(female) + as.factor(cohort_group) + baseline_wealth + baseline_income + edu_years + center_age * race, data=model_hrs, weights=pweight)
library(splines)
model1 <- lmerTest::lmer(formula1, weights=pweight, data=model_hrs, REML=FALSE)
model2 <- lmerTest::lmer(formula2, weights=pweight, data=model_hrs, REML=FALSE)
model3 <- lmerTest::lmer(formula3, weights=pweight, data=model_hrs, REML=FALSE)
# model3_spline <- lmerTest::lmer(formula3_spline, weights=pweight, data=model_hrs, REML=FALSE)
# model4 <- lmerTest::lmer(formula4, weights=pweight, data=model_hrs, REML=FALSE)
# model5 <- lmerTest::lmer(formula5, weights=pweight, data=model_hrs, REML=FALSE)

test <- data.table(y=c(1:25,25:1), age=1:50)
m <- lm(y ~ bs(age, degree=1, knots=c(25)), data=test)
test[, pred := predict(m,data=test)]
ggplot(data=test) +
  geom_point(aes(x=age,y=y)) +
  geom_line(aes(x=age,y=pred))
24.69065/(25 - 1)
(0.47002-24.69065)/(50 - 25)

## CONVERT SPLINE ESTIMATES TO SLOPE COEFFICIENTS
## https://stackoverflow.com/questions/37362738/how-to-interpret-lm-coefficient-estimates-when-using-bs-function-for-splines
# knot.boundary.left <- min(model_hrs[, center_age])
# knot.boundary.right <- max(model_hrs[, center_age])
# knot <- 4
# fit <- model3
# slope.1 <- summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)1',1] /(knot - knot.boundary.left)
# slope.2 <- (summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)2',1] - summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)1',1]) / (knot.boundary.right - knot)

## TRY SPLINE PREDICTIONS
# spline_table[, center_age := center_age * 5 + 50]
# spline_table[, female := 0]
# spline_table[, cohort_group := '1930']
# spline_table[, edu_years := 0]
# spline_table[, ('baseline_cog') := 0]
# spline_table[, ('baseline_wealth') := 0]
# spline_table[, ('baseline_income') := 0]
# library(merTools)
# spline_table[, id_factor := '502223020']
# spline_table[, matchvec_1 := 1]
# test <- predictInterval(model3_spline, which='fixed', level=0.95, newdata = spline_table, n.sims = 1)
# spline_table[, pred_mean := predict(model3_spline, re.form=NA, spline_table)]
# ggplot(data=spline_table, aes(x=center_age,y=pred_mean)) + 
#   # geom_ribbon(aes(ymin=pred_lower,
#   #                 ymax=pred_upper,
#   #                 fill=race),
#   #             color='black',
#   #             alpha=0.6,
#   #             size=0.2) + 
#   geom_line(aes(group=race),color='black', size=1) + 
#   geom_vline(xintercept = 69.9, size=2) + 
#   theme_bw() + 
#   labs(y='Predicted mean cognition score (standard deviations)',x='Age',fill='Race') + 
#   theme(strip.text.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20, margin = margin(r=10)),
#         axis.title.x = element_text(size = 20, margin = margin(t=10)),
#         axis.text = element_text(size = 20),
#         legend.key.size = unit(3,'line'),
#         legend.title = element_text(size = 20),
#         legend.text = element_text(size = 15))

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
# d2 <- merge(model_hrs, d[!is.na(matchvec_1), c('id_factor','matchvec_1')], by='id_factor')
d2 <- merge(model_hrs, d[, c('id_factor','matchvec_1')], by='id_factor', all.x = T)
hrs_no_match <- copy(d2[race=='white',])
hrs_no_match[, race := ifelse(is.na(matchvec_1), 'white', 'black')]
d2 <- d2[!is.na(matchvec_1),]
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
pdf(paste0('matched_pairs_edu_cat', file_tag, '.pdf'), width=12, height=8)
print(gg.love)
print(gg.race)
print(gg.race.edu)
dev.off()

## Try modelling with fixed effect on matched pairs.
# model6 <- lmerTest::lmer(formula5, weights=pweight, data=d2, REML=FALSE)
formula4 <- as.formula(paste0(outcome_var, ' ~ as.factor(matchvec_1) + center_age*over_70*race + center_age*edu_years + (center_age||id_factor)'))
formula4_spline <- as.formula(paste0(outcome_var, ' ~ as.factor(matchvec_1) + bs(center_age, knots=4, degree=1)*race +
  bs(center_age, knots=4, degree=1)*edu_years + (center_age||id_factor)'))

d2[, race := factor(race, levels=c('white','black'))]
model4 <- lmerTest::lmer(formula4, weights=pweight, data=d2, REML=FALSE)
model4_spline <- lmerTest::lmer(formula4_spline, weights=pweight, data=d2[!is.na(matchvec_1)], REML=FALSE)

## Make model table
# tab_model(list(model2,model3), dv.labels = c('Model 1','Model 2','Matched'), show.ci=FALSE, p.style='asterisk', file='three_way_interaction_edu_years.html')

## Try making predicted scores
coef_table <- as.data.table(expand.grid(seq(0,7,0.1), c('white','black')))
setnames(coef_table, c('center_age','race'))
coef_table[, raceblack := ifelse(race=='black',1,0)]
coef_table[, over_70 := ifelse(center_age>=4, 1, 0)]
## TARGET AVERAGE: MALE, 1930 COHORT, AVERAGE WEALTH/INCOME/EDUCATION
## EDU * AGE DUMMIES
  # coef_table[, (paste0('center_age:edu_years')) := center_age * edu_years]
## over_70 * RACE DUMMIES
  coef_table[, (paste0('over_70:raceblack')) := over_70 * raceblack]
## over_70 * AGE DUMMIES
  coef_table[, ('center_age:over_70') := over_70 * center_age]
## RACE * AGE DUMMIES
  coef_table[, ('center_age:raceblack') := center_age * raceblack]
## over_70 * AGE * RACE DUMMIES
  coef_table[, (paste0('center_age:over_70:raceblack')) := center_age * over_70 * raceblack]
## Other dummies in case
  coef_table[, ('as.factor(female)1') := 0]
  if(oldest) {
    coef_table[, ('as.factor(cohort_group)1930') := 0]
  }
  if(!(oldest)) {
  coef_table[, ('as.factor(cohort_group)1930') := 0]
  coef_table[, ('as.factor(cohort_group)1940') := 0]
  coef_table[, ('as.factor(cohort_group)1950') := 0]
  }
  coef_table[, ('baseline_cog') := 0]
  coef_table[, ('baseline_wealth') := 0]
  coef_table[, ('baseline_income') := 0]
coef_table[, ("(Intercept)") := 1]
spline_table <- copy(coef_table)
betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway.
betas <- betas[, !c('matchvec_1','edu_years','center_age:edu_years'), with=F]
beta_names <- names(betas)
for(d in 1:1000) {
  betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
  names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway.
  betas <- betas[, !c('matchvec_1','edu_years','center_age:edu_years'), with=F]
  beta_names <- names(betas)
  template <- coef_table[, beta_names, with=F]
  setcolorder(template, beta_names)
  coef_table[, (paste0('draw',d)) := as.numeric(as.matrix(betas) %*% t(as.matrix(template)))]
}
coef_table[, pred_mean := apply(.SD,1,mean), .SDcols=paste0('draw',1:1000)]
coef_table[, pred_lower := apply(.SD,1,quantile,probs=0.025), .SDcols=paste0('draw',1:1000)]
coef_table[, pred_upper := apply(.SD,1,quantile,probs=0.975), .SDcols=paste0('draw',1:1000)]
# coef_table[, race_edu := paste0(race,'_',edu_cat)]
coef_table[, center_age := center_age * 5 + 50]
coef_table[race=='black', race := 'Black']
coef_table[race=='white', race := 'White']
# coef_table[edu_years==-1.903866, edu_years_string := '8 years']
# coef_table[edu_years==-0.3609927, edu_years_string := '12 years']
# coef_table[edu_years==1.181881, edu_years_string := '16 years']
# coef_table[, edu_years_string := factor(edu_years_string, levels=c('8 years','12 years','16 years'))]
pdf(paste0('C:/Users/ngraetz/Documents/repos/hrs/Figure1_Model3_', file_tag, '_ageinteraction.pdf'), width=12, height=8)
ggplot(data=coef_table, aes(x=center_age,y=pred_mean)) + 
  geom_ribbon(aes(ymin=pred_lower,
                  ymax=pred_upper,
                  fill=race),
              color='black',
              alpha=0.6,
              size=0.2) + 
  geom_line(aes(group=race),color='black', size=1) + 
  geom_vline(xintercept = 69.9, size=2) + 
  # facet_wrap(~edu_years_string) + 
  theme_bw() + 
  labs(y='Predicted mean cognition score (standard deviations)',x='Age',fill='Race') + 
  theme(strip.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.text = element_text(size = 20),
        legend.key.size = unit(3,'line'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
dev.off()

## TRY SPLINE PREDICTIONS
spline_table[, center_age := center_age * 5 + 50]
spline_table[, female := 0]
spline_table[, cohort_group := '1930']
spline_table[, edu_years := 0]
spline_table[, ('baseline_cog') := 0]
spline_table[, ('baseline_wealth') := 0]
spline_table[, ('baseline_income') := 0]
library(merTools)
spline_table[, id_factor := '502223020']
spline_table[, matchvec_1 := 1]

test <- predictInterval(model3_spline, which='full', level=0.95, newdata = spline_table, n.sims = 1)
spline_table[, pred_mean := predict(model3_spline, re.form=NA, spline_table)]

test <- predict(model4_spline, spline_table)
spline_table[, pred_mean := test]

# test <- predictInterval(model3_spline, which='full', level=0.95, newdata = spline_table, n.sims = 1)
# spline_table[, pred_mean := predict(model3_spline, re.form=NA, spline_table)]
pdf(paste0('C:/Users/ngraetz/Documents/repos/hrs/Figure1_Model3_', file_tag, '_ageinteraction_spline.pdf'), width=12, height=8)
ggplot(data=spline_table, aes(x=center_age,y=pred_mean)) + 
  # geom_ribbon(aes(ymin=pred_lower,
  #                 ymax=pred_upper,
  #                 fill=race),
  #             color='black',
  #             alpha=0.6,
  #             size=0.2) + 
  geom_line(aes(color=race), size=2) + 
  geom_vline(xintercept = 70, size=1) + 
  # facet_wrap(~edu_years_string) + 
  theme_bw() + 
  labs(y='Predicted mean cognition score (standard deviations)',x='Age',fill='Race') + 
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
pdf(paste0('C:/Users/ngraetz/Documents/repos/hrs/predicted_coefficients_', file_tag, '.pdf'), width=8, heigh=6)
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

## NEW TABLE 1
library(survey)
table1 <- copy(model_hrs)
table1 <- copy(hrs_no_match)
table1[, original_baseline_income := exp(original_baseline_income)]
table1[, male := ifelse(female==0, 1, 0)]
table1[, edu_years := as.numeric(educYrs)]
for(var in c('cohort_group','age_group')) {
  for(v in unique(table1[, get(var)])) {
    table1[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}
table1[, race := as.character(race)]
des <- svydesign(id=~id, weight=~pweight, data=table1[!is.na(cognitive), ])
get_table1 <- function(v) {
  p <- svyttest(formula=as.formula(paste0(v,'~race')), design=des)$p.value
  p <- data.table(type='pvalue',mean=p,se=0)
  race_means <- as.data.table(svyby(as.formula(paste0('~',v)), ~race, design=des, svymean, na.rm=T))
  setnames(race_means, c('type','mean','se'))
  full_mean <- as.data.table(svymean(as.formula(paste0('~',v)), design=des, na.rm=T))
  full_mean[, type := 'full']
  setnames(full_mean, c('mean','se','type'))
  all <- rbind(race_means, full_mean, p)
  all[, variable := v]
  return(all)
}
table_vars <- c('original_cognitive','original_baseline_cog','income','original_baseline_income','wealth',
                'original_baseline_wealth','edu_years','male','female','cohort_group_1930','cohort_group_1940',
                'cohort_group_1950','age')
for(v in table_vars) table1 <- table1[!is.na(get(v)), ]
table1 <- rbindlist(lapply(table_vars, get_table1))
table1[, mean := round(mean, 2)]
table1[, se := round(se, 2)]
table1[variable=='income', mean := (mean/1000)]
table1[variable=='income', se := (se/1000)]
table1[variable=='wealth', mean := (mean/1000)]
table1[variable=='wealth', se := (se/1000)]
table1[variable=='original_baseline_wealth', mean := (mean/1000)]
table1[variable=='original_baseline_wealth', se := (se/1000)]
table1[variable=='original_baseline_income', mean := (mean/1000)]
table1[variable=='original_baseline_income', se := (se/1000)]
table1 <- dcast(table1, variable ~ type, value.var = c('mean','se'))
table1 <- table1[, c('variable','mean_full','se_full','mean_white','se_white','mean_black','se_black','mean_pvalue'), with=F]
setcolorder(table1, c('variable','mean_full','se_full','mean_white','se_white','mean_black','se_black','mean_pvalue'))
table1[, p := '']
table1[mean_pvalue<=0.05, p := paste0(p,'*')]
table1[mean_pvalue<=0.01, p := paste0(p,'*')]
table1[mean_pvalue<=0.001, p := paste0(p,'*')]
table1[, mean_pvalue := NULL]
for(v in names(table1)[!(names(table1) %in% c('variable','p'))]) table1[, (v) := format(get(v), nsmall=0, digits=2, big.mark=",")]
for(v in c('se_full','se_white','se_black')) table1[, (v) := paste0('(',get(v),')')]
clean_names <- data.table(variable=c('original_cognitive','original_baseline_cog','income',
                                     'original_baseline_income','wealth','original_baseline_wealth',
                                     'edu_years','male','female','cohort_group_1930','cohort_group_1940',
                                     'cohort_group_1950','age',paste0('age_group_',seq(50,75,5))),
                          name=c('Cognition (all ages)','Cognition (baseline)','Income (all ages)','Income (baseline)',
                                 'Wealth (all ages)','Wealth (baseline)','Years of education','Male','Female','1930 cohort',
                                 '1940 cohort','1950 cohort','Age','Age 50-54','Age 55-59','Age 60-64','Age 65-69','Age 70-74',
                                 'Age 75-79'))
clean_names[, cov_sort := 1:.N]
table1 <- merge(table1, clean_names, by='variable')
table1 <- table1[order(cov_sort)]
table1[, cov_sort := NULL]
table1[, variable := NULL]
setcolorder(table1, 'name')

first_header <- as.list(c('',rep(c('Mean','SE'),3),''))
second_header <- as.list(c('',rep('Full Sample',2),rep('White Sample',2),rep('Black Sample',2),'Difference'))
names(first_header) <- names(table1)
names(second_header) <- names(table1)

library(flextable)
library(officer)
ft <- flextable(table1, theme_fun = theme_booktabs) %>%
  delete_part(part = "header") %>%
  add_header(values=first_header) %>%
  add_header(values=second_header) %>%
  merge_h(part = "header") %>%
  align(align = 'right', part='body') %>%
  align(j=1, align = 'left', part='body') %>%
  align(i=1, align = 'right', part='header') %>%
  align(i=2, align = 'right', part='header') %>%
  hline_top(part = 'body', border = fp_border(color='black',width=2)) %>%
  hline_top(part = 'header', border = fp_border(color='black',width=2)) %>%
  hline(j=2:8,part = 'header', border = fp_border(color='black',width=2)) %>%
  padding(padding = 0.0001, part = 'all') %>%
  padding(padding.top = 0, part='all') %>%
  # width(j=2:8, width=0.4) %>%
  # width(j=1, width=2.23) %>%
  add_footer_lines(c('* p < 0.05; ** p < 0.01; *** p < 0.001')) %>%
  fontsize(size = 12, part='all') %>%
  font(fontname = 'Times New Roman', part='all') %>%
  autofit()

setwd('C:/Users/ngraetz/Documents/repos/hrs/')
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print(target = paste0("table1", file_tag, "_matched_comparison.docx"))

## Table 2: models
## Use summary() in lmerTest for p-values and r.squaredGLMM() in MuMIn for R^2
library(lmerTest)
library(MuMIn)
get_coefs <- function(x) {
  this_sum <- summary(get(x))
  coefs <- as.data.table(this_sum$coefficients)
  coefs[, variable := rownames(this_sum$coefficients)]
  coefs <- coefs[, c('Estimate','Std. Error','Pr(>|t|)','variable')]
  setnames(coefs, c('mean','se','pvalue','variable'))
  coefs[, p := '']
  coefs[pvalue<=0.05, p := paste0(p,'*')]
  coefs[pvalue<=0.01, p := paste0(p,'*')]
  coefs[pvalue<=0.001, p := paste0(p,'*')]
  coefs[, pvalue := NULL]
  rs <- r.squaredGLMM(get(x))
  rs <- data.table(variable=c('mar_rsq','cond_rsq'),mean=c(rs[1],rs[2]))
  bic <- data.table(variable=c('bic'),mean=round(this_sum$AICtab['BIC']))
  py <- data.table(variable=c('N'),mean=this_sum$devcomp$dims['N'])
  people <- data.table(variable=c('people'),mean=this_sum$ngrps)
  coefs <- rbindlist(list(coefs, rs, bic, py, people), fill=T)
  coefs[, model := x]
  coefs[variable=="center_age:edu_years", variable := "edu_years:center_age"]
  coefs[variable=="center_age:edu_years:raceblack", variable := "edu_years:center_age:raceblack"]
  coefs[variable=='center_age:raceblack:edu_years_mother', variable := 'center_age:edu_years_mother:raceblack']
  coefs[variable=="raceblack:edu_years_mother", variable := 'edu_years_mother:raceblack"']
  return(coefs)
}
options(scipen = 999)
coefs <- rbindlist(lapply(c('model1','model2','model3','model4'), get_coefs))
coefs[is.na(p), p := '']
for(v in c('mean','se')) coefs[, (v) := format(round(get(v),2), nsmall=0, big.mark=",")]
coefs[se=='     NA', se := '']
clean_names <- data.table(variable=c("(Intercept)","as.factor(female)1","as.factor(cohort_group)1940","as.factor(cohort_group)1950","baseline_wealth","baseline_income","baseline_cog","edu_years","edu_years_mother","center_age","raceblack",
                                     "edu_years:center_age",
                                     "center_age:raceblack",
                                     "edu_years:raceblack",
                                     "edu_years:center_age:raceblack",
                                     "center_age:edu_years_mother",
                                     "edu_years_mother:raceblack",
                                     "center_age:edu_years_mother:raceblack",
                                     'over_70','over_70:raceblack','center_age:over_70',
                                     "center_age:over_70:raceblack",
                                     "mar_rsq","cond_rsq","bic","N","people"),
                          name=c('Intercept','Female','1940 cohort','1950 cohort','Wealth (baseline)','Income (baseline)','Cognition (baseline)','Years of education','Mother years of education','Age','Black','Age*Edu','Age*Black','Edu*Black','Age*Edu*Black','Age*MEdu','MEdu*Black','Age*MEdu*Black',
                                 'Over70','Over70*Black','Age*Over70','Age*Over70*Black',
                                 'Marginal R2','Conditional R2','BIC','Person-years','Individuals'))
clean_names[, cov_sort := 1:.N]
coefs <- merge(coefs, clean_names, by='variable', all.x=T)
coefs <- coefs[!is.na(name), ]
coefs <- dcast(coefs, cov_sort+name~model, value.var = c('mean','se','p'))
coefs <- coefs[order(cov_sort)]
coefs[, c('cov_sort') := NULL]
for(v in names(coefs)[names(coefs)!='name']) coefs[is.na(get(v)), (v) := '']
for(v in paste0('se_model',c(1:4))) {
  coefs[get(v)=='  NA', (v) := '']
  coefs[get(v)!='', (v) := paste0('(',get(v),')')]
}
setcolorder(coefs, c('name','mean_model1','se_model1','p_model1',
                     'mean_model2','se_model2','p_model2',
                     'mean_model3','se_model3','p_model3',
                     'mean_model4','se_model4','p_model4'))

first_header <- as.list(c('',rep(c('Coef.','SE',''),4)))
second_header <- as.list(c('',rep('Model 1',3),rep('Model 2',3),rep('Model 3',3),rep('Model 4',3)))
names(first_header) <- names(coefs)
names(second_header) <- names(coefs)

library(flextable)
library(officer)
ft <- flextable(coefs, theme_fun = theme_booktabs) %>%
  delete_part(part = "header") %>%
  add_header(values=first_header) %>%
  add_header(values=second_header) %>%
  merge_h(part = "header") %>%
  align(align = 'right', part='body') %>%
  align(j=1, align = 'left', part='body') %>%
  align(i=1, align = 'right', part='header') %>%
  align(i=2, align = 'right', part='header') %>%
  hline_top(part = 'body', border = fp_border(color='black',width=2)) %>%
  hline_top(part = 'header', border = fp_border(color='black',width=2)) %>%
  hline(j=2:13,part = 'header', border = fp_border(color='black',width=2)) %>%
  hline(j=1:13,i=16,part = 'body', border = fp_border(color='black',width=2)) %>%
  padding(padding = 0.0001, part = 'all') %>%
  padding(padding.top = 0, part='all') %>%
  # width(j=2:8, width=0.4) %>%
  # width(j=1, width=2.23) %>%
  add_footer_lines(c('* p < 0.05; ** p < 0.01; *** p < 0.001')) %>%
  fontsize(size = 12, part='all') %>%
  font(fontname = 'Times New Roman', part='all') %>%
  autofit()

setwd('C:/Users/ngraetz/Documents/repos/hrs/')
doc <- read_docx() %>%
  body_add_flextable(value = ft, split = TRUE) %>%
  body_end_section_landscape() %>% # a landscape section is ending here
  print(target = paste0("table2_", file_tag, "_over70interaction.docx"))
