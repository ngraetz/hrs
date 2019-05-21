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

outcome_var <- 'cognitive'
REML <- FALSE
weight <- FALSE
ses <- FALSE
if(ses) numeric_vars <- c('edu_years','wealth','log_income') 
if(!ses) numeric_vars <- NULL

numeric_vars <- 'edu_years'
factor_vars <- c('as.factor(female)','as.factor(cohort_group)')

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

## Subset to respondents with at least 3 responses on outcome.
id_nums <- model_hrs[, .N , by = id]
dim(id_nums[N == 1, ])
model_hrs <- model_hrs[id %in% id_nums[N > 3, id], ]
# model_hrs[, id_factor := as.factor(id)]

## Rescale all variables to mean/sd=0/1. Keep track of means/sds to transform back after predicting.
scale_vars <- c(outcome_var, numeric_vars, 'baseline_cog')
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

formula1 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars), collapse=' + '), ' + center_age*baseline_cog*race + (center_age||id_factor)'))
model1 <- lmer(formula1, weights=pweight, data=model_hrs, REML=FALSE)

formula2 <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars), collapse=' + '), ' + baseline_cog + center_age*edu_years*race + (center_age||id_factor)'))
model2 <- lmer(formula2, weights=pweight, data=model_hrs, REML=FALSE)

## TRY MATCHING BASELINE OBSERVATIONS 
source('C:/Users/ngraetz/Documents/repos/hrs/causal_functions.R')
d <- model_hrs[age==age_obs & race %in% c('black','white') & !is.na(edu_years), ]
for(var in c('cohort_group','race')) {
  for(v in unique(d[, get(var)])) {
    d[, (paste0(var,'_',v)) := ifelse(get(var)==v, 1, 0)]
  }
}
confounders <- c('age', 'female', 'cohort_group_1930', 'cohort_group_1940', 'baseline_cog', 'edu_years', 'wealth', 'log_income')
treatment <- 'race_black'

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
pdf('matched_pairs.pdf', width=12, height=8)
print(gg.love)
print(gg.race)
print(gg.race.edu)
dev.off()

## Try modelling with fixed effect on matched pairs.
formula3 <- as.formula(paste0(outcome_var, ' ~ matchvec_1 + center_age*edu_years*race + (center_age||id_factor)'))
model3 <- lmer(formula3, weights=pweight, data=d2, REML=FALSE)

## Make model table
tab_model(list(model1,model2,model3), dv.labels = c('Model 1','Model 2','Matched'), show.ci=FALSE, p.style='asterisk', file='three_way_interaction.html')

