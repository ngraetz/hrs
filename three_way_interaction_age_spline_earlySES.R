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
library(splines)

repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
hrs <- readRDS('hrs_imputed_v3.RDS')
oldest <- FALSE
drop_oldest <- FALSE
file_tag <- 'age_interaction_earlySES'

## Get early life
early <- fread("C:/Users/ngraetz/Downloads/rand_childhood_hrs.csv")

## Restrict race
hrs <- hrs[race %in% c('white','black'), ]
hrs[, race := factor(race, levels=c('white','black'))]
hrs[, edu_cat := factor(edu_cat, levels=c('highschool','less_highschool','some_college','college'))]
hrs <- merge(hrs, early, by.x='id', by.y='rahhidpn')
hrs[, mother_edu := rameduc]
hrs[, child_health := chealth]

outcome_var <- 'cognitive'
REML <- FALSE
weight <- FALSE

# 'baseline_wealth','baseline_income','baseline_cognition'
numeric_vars <- c('edu_years','mother_edu','child_health','baseline_wealth','baseline_income','baseline_cog') 
factor_vars <- c('as.factor(female)','as.factor(cohort_group)')

## MATCHING
ifelse(oldest, confounders <- c('age', 'female', 'cohort_group_1930', 'baseline_cog', 'edu_years', 'baseline_wealth', 'baseline_income'), confounders <- c('age', 'female', 'cohort_group_1930', 'cohort_group_1940', 'baseline_cog', 'edu_years','mother_edu','child_health','baseline_wealth', 'baseline_income'))
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
scale_vars <- c(outcome_var, numeric_vars) ## Drop outcome from scale vars for now
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

for(v in c('cognitive','baseline_cog','baseline_wealth','baseline_income','edu_years_mother','edu_years','mother_edu','child_health','cohort_group')) model_hrs <- model_hrs[!is.na(get(v)), ]
model_hrs[, demean_cognitive := mean(cognitive), by='id']
model_hrs[, demean_cognitive := cognitive - demean_cognitive]

formula1 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
  bs(center_age, knots=4, degree=1)*race +
  (center_age||id_factor)
formula2 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
  bs(center_age, knots=4, degree=1)*race +
  baseline_cog + 
  (center_age||id_factor)
# formula3 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   center_age*edu_years_mother*race + 
#   baseline_cog + (center_age||id_factor)
formula3 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
  bs(center_age, knots=4, degree=1)*race +
  bs(center_age, knots=4, degree=1)*edu_years + 
  baseline_cog + baseline_wealth + baseline_income + mother_edu + child_health + 
  (center_age||id_factor)
# formula3_spline <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   bs(center_age, knots=4, degree=1)*race +
#   bs(center_age, knots=4, degree=1)*edu_years + 
#   baseline_cog + baseline_wealth + baseline_income + mother_edu + child_health + 
#   (center_age||id_factor)
# formula5 <- cognitive ~ as.factor(female) + as.factor(cohort_group) + 
#   center_age*edu_years*race + center_age*edu_years_mother*race + 
#   baseline_cog + baseline_wealth + baseline_income + (center_age||id_factor)

# ols1 <- lm(demean_cognitive ~ as.factor(female) + as.factor(cohort_group) + baseline_wealth + baseline_income + edu_years + center_age * race, data=model_hrs, weights=pweight)
model1 <- lmerTest::lmer(formula1, weights=pweight, data=model_hrs, REML=FALSE)
model2 <- lmerTest::lmer(formula2, weights=pweight, data=model_hrs, REML=TRUE)
model3 <- lmerTest::lmer(formula3, weights=pweight, data=model_hrs, REML=FALSE)
# model3_spline <- lmerTest::lmer(formula3_spline, weights=pweight, data=model_hrs, REML=FALSE)
# model4 <- lmerTest::lmer(formula4, weights=pweight, data=model_hrs, REML=FALSE)
# model5 <- lmerTest::lmer(formula5, weights=pweight, data=model_hrs, REML=FALSE)

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
formula4 <- as.formula(paste0(outcome_var, ' ~ as.factor(matchvec_1) + bs(center_age, knots=4, degree=1)*race + bs(center_age, knots=4, degree=1)*edu_years + (center_age||id_factor)'))
# formula4_spline <- as.formula(paste0(outcome_var, ' ~ as.factor(matchvec_1) + bs(center_age, knots=4, degree=1)*race + bs(center_age, knots=4, degree=1)*edu_years + (center_age||id_factor)'))
d2[, race := factor(race, levels=c('white','black'))]
model4 <- lmerTest::lmer(formula4, weights=pweight, data=d2, REML=FALSE)
# model4_spline <- lmerTest::lmer(formula4_spline, weights=pweight, data=d2, REML=FALSE)

# knot.boundary.left <- min(d2[, center_age])
# knot.boundary.right <- max(d2[, center_age])
# knot <- 4
# fit <- model4_spline
# slope.1 <- summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)1',1] /(knot - knot.boundary.left)
# slope.2 <- (summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)2',1] - summary(fit)$coefficients['bs(center_age, knots = 4, degree = 1)1',1]) / (knot.boundary.right - knot)

## Try making predicted scores
coef_table <- as.data.table(expand.grid(seq(0,5.2,0.1), c('white','black')))
setnames(coef_table, c('center_age','race'))
coef_table[, raceblack := ifelse(race=='black',1,0)]
# coef_table[, over_70 := ifelse(center_age>=4, 1, 0)]
coef_table[, ('bs(center_age, knots = 4, degree = 1)1') := bs(center_age,knots=4,degree=1)[,1]]
coef_table[, ('bs(center_age, knots = 4, degree = 1)2') := bs(center_age,knots=4,degree=1)[,2]]
## TARGET AVERAGE: MALE, 1930 COHORT, AVERAGE WEALTH/INCOME/EDUCATION
## EDU * AGE DUMMIES (AT AVERAGE, EDU=0)
# coef_table[, (paste0('bs(center_age, knots = 4, degree = 1)1:edu_years')) := get('bs(center_age, knots = 4, degree = 1)1') * edu_years]
# coef_table[, (paste0('bs(center_age, knots = 4, degree = 1)2:edu_years')) := get('bs(center_age, knots = 4, degree = 1)2') * edu_years]
# ## over_70 * RACE DUMMIES
# coef_table[, (paste0('over_70:raceblack')) := over_70 * raceblack]
# ## over_70 * AGE DUMMIES
# coef_table[, ('center_age:over_70') := over_70 * center_age]
## RACE * AGE DUMMIES
coef_table[, ('bs(center_age, knots = 4, degree = 1)1:raceblack') := get('bs(center_age, knots = 4, degree = 1)1')*raceblack]
coef_table[, ('bs(center_age, knots = 4, degree = 1)2:raceblack') := get('bs(center_age, knots = 4, degree = 1)2')*raceblack]
## over_70 * AGE * RACE DUMMIES
# coef_table[, (paste0('center_age:over_70:raceblack')) := center_age * over_70 * raceblack]
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
coef_table[, ('mother_edu') := 0]
coef_table[, ('child_health') := 0]
coef_table[, ("(Intercept)") := 1]
spline_table <- copy(coef_table)
betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway
betas <- betas[, !c(names(betas)[!(names(betas) %in% names(coef_table))]), with=F]
beta_names <- names(betas)
all_betas <- arm::sim(model3, n.sim=1000)

for(d in 1:1000) {
  betas <- as.data.table(all_betas@'fixef')[d,]
  betas <- betas[, !c(names(betas)[!(names(betas) %in% names(coef_table))]), with=F]
  beta_names <- names(betas)
  template <- coef_table[, beta_names, with=F]
  setcolorder(template, beta_names)
  coef_table[, (paste0('draw',d)) := as.numeric(as.matrix(betas) %*% t(as.matrix(template)))]
  ## Convert to actual score space
  coef_table[, (paste0('draw',d)) := (get(paste0('draw',d))*scales[variable=='cognitive', sd]) + scales[variable=='cognitive', mean]]
}
coef_table[, pred_mean := apply(.SD,1,mean), .SDcols=paste0('draw',1:400)]
coef_table[, pred_lower := apply(.SD,1,quantile,probs=0.025), .SDcols=paste0('draw',1:400)]
coef_table[, pred_upper := apply(.SD,1,quantile,probs=0.975), .SDcols=paste0('draw',1:400)]
coef_table[, female := as.factor(0)]
coef_table[, cohort_group := as.factor(1930)]
coef_table[, edu_years := 0]
coef_table[, id_factor := '502223020']
coef_table[, pred_full := predict(model3, re.form=~0, newdata=coef_table)]
coef_table[, pred_full := (pred_full*scales[variable=='cognitive', sd]) + scales[variable=='cognitive', mean]]
coef_table[, center_age := center_age * 5 + 50]
coef_table[race=='black', race := 'Black']
coef_table[race=='white', race := 'White']
# pdf(paste0('C:/Users/ngraetz/Documents/repos/hrs/Figure1_Model3_', file_tag, '.pdf'), width=12, height=8)
gg_pred_race <- ggplot(data=coef_table, aes(x=center_age,y=pred_mean)) + 
  geom_ribbon(aes(ymin=pred_lower,
                  ymax=pred_upper,
                  fill=race),
              color='black',
              alpha=0.6,
              size=0.2) + 
  geom_line(data=coef_table,aes(x=center_age,y=pred_full,group=race),color='black',size=1) + 
  geom_vline(xintercept = 70, size=1) + 
  # facet_wrap(~edu_years_string) + 
  theme_bw() + 
  labs(y='Predicted mean cognition score',x='Age',fill='Race') + 
  theme(strip.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.text = element_text(size = 20),
        legend.key.size = unit(3,'line'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
# dev.off()

## Try making predicted scores
model_hrs[, edu_actual := as.numeric(educYrs)]
model_hrs[, edu_back := (edu_years*scales[variable=='edu_years', sd])+scales[variable=='edu_years', mean]]
edu_levels <- unique(model_hrs[edu_back %in% c(8,12,16), c('edu_years')])
coef_table <- as.data.table(expand.grid(seq(0,5.2,0.1), c(edu_levels$edu_years)))
setnames(coef_table, c('center_age','edu_years'))
coef_table[, race := 'black']
coef_table[, raceblack := ifelse(race=='black',1,0)]
# coef_table[, over_70 := ifelse(center_age>=4, 1, 0)]
coef_table[, ('bs(center_age, knots = 4, degree = 1)1') := bs(center_age,knots=4,degree=1)[,1]]
coef_table[, ('bs(center_age, knots = 4, degree = 1)2') := bs(center_age,knots=4,degree=1)[,2]]
## TARGET AVERAGE: MALE, 1930 COHORT, AVERAGE WEALTH/INCOME/EDUCATION
## EDU * AGE DUMMIES 
coef_table[, (paste0('bs(center_age, knots = 4, degree = 1)1:edu_years')) := get('bs(center_age, knots = 4, degree = 1)1') * edu_years]
coef_table[, (paste0('bs(center_age, knots = 4, degree = 1)2:edu_years')) := get('bs(center_age, knots = 4, degree = 1)2') * edu_years]
# ## over_70 * RACE DUMMIES
# coef_table[, (paste0('over_70:raceblack')) := over_70 * raceblack]
# ## over_70 * AGE DUMMIES
# coef_table[, ('center_age:over_70') := over_70 * center_age]
## RACE * AGE DUMMIES
coef_table[, ('bs(center_age, knots = 4, degree = 1)1:raceblack') := get('bs(center_age, knots = 4, degree = 1)1')*raceblack]
coef_table[, ('bs(center_age, knots = 4, degree = 1)2:raceblack') := get('bs(center_age, knots = 4, degree = 1)2')*raceblack]
## over_70 * AGE * RACE DUMMIES
# coef_table[, (paste0('center_age:over_70:raceblack')) := center_age * over_70 * raceblack]
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
coef_table[, ('mother_edu') := 0]
coef_table[, ('child_health') := 0]
coef_table[, ("(Intercept)") := 1]
spline_table <- copy(coef_table)
betas <- as.data.table(t(MASS::mvrnorm(1, mu = fixef(model3), Sigma = vcov(model3))))
names(betas)[!(names(betas) %in% names(coef_table))] ## Don't predict with "matchvec_1", should be virtually 0 anyway
betas <- betas[, !c(names(betas)[!(names(betas) %in% names(coef_table))]), with=F]
beta_names <- names(betas)
all_betas <- arm::sim(model3, n.sim=1000)

for(d in 1:1000) {
  betas <- as.data.table(all_betas@'fixef')[d,]
  # betas <- betas[, !c(names(betas)[!(names(betas) %in% names(coef_table))]), with=F]
  beta_names <- names(betas)
  template <- coef_table[, beta_names, with=F]
  setcolorder(template, beta_names)
  coef_table[, (paste0('draw',d)) := as.numeric(as.matrix(betas) %*% t(as.matrix(template)))]
  ## Convert to actual score space
  coef_table[, (paste0('draw',d)) := (get(paste0('draw',d))*scales[variable=='cognitive', sd]) + scales[variable=='cognitive', mean]]
}
coef_table[, pred_mean := apply(.SD,1,mean), .SDcols=paste0('draw',1:400)]
coef_table[, pred_lower := apply(.SD,1,quantile,probs=0.025), .SDcols=paste0('draw',1:400)]
coef_table[, pred_upper := apply(.SD,1,quantile,probs=0.975), .SDcols=paste0('draw',1:400)]
coef_table[, female := as.factor(0)]
coef_table[, cohort_group := as.factor(1930)]
coef_table[, id_factor := '502223020']
coef_table[, pred_full := predict(model3, re.form=~0, newdata=coef_table)]
coef_table[, pred_full := (pred_full*scales[variable=='cognitive', sd]) + scales[variable=='cognitive', mean]]
coef_table[, center_age := center_age * 5 + 50]
coef_table[race=='black', race := 'Black']
coef_table[race=='white', race := 'White']
coef_table[edu_years < -1, edu_level := '8 years']
coef_table[edu_years > 1, edu_level := '16 years']
coef_table[is.na(edu_level), edu_level := '12 years']
coef_table[, edu_level := factor(edu_level, levels=c('8 years','12 years','16 years'))]
gg_pred_edu <- ggplot(data=coef_table, aes(x=center_age,y=pred_mean)) + 
  geom_ribbon(aes(ymin=pred_lower,
                  ymax=pred_upper,
                  fill=as.factor(edu_level)),
              color='black',
              alpha=0.6,
              size=0.2) + 
  geom_line(data=coef_table,aes(x=center_age,y=pred_full,group=as.factor(edu_level)),color='black',size=1) + 
  geom_vline(xintercept = 70, size=1) + 
  # facet_wrap(~edu_years_string) + 
  theme_bw() + 
  labs(y='Predicted mean cognition score',x='Age',fill='Education') + 
  theme(strip.text.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r=10)),
        axis.title.x = element_text(size = 20, margin = margin(t=10)),
        axis.text = element_text(size = 20),
        legend.key.size = unit(3,'line'),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
pdf(paste0('C:/Users/ngraetz/Documents/repos/hrs/Figure1_predicted_socres_', file_tag, '.pdf'), width=12, height=16)
grid.arrange(gg_pred_race,gg_pred_edu,ncol=1)
dev.off()

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
  ## Drop pair id fixed effects from matched model
  coefs <- coefs[!grepl('matchvec_1',variable),]
  ## Convert spline terms (mean/SE) to linear slopes
  ## https://stackoverflow.com/questions/37362738/how-to-interpret-lm-coefficient-estimates-when-using-bs-function-for-splines/37363943
  knot.boundary.left <- 0
  knot.boundary.right <- 5.2
  knot <- 4
  draws <- arm::sim(get(x),n.sim=1000)
  draws <- as.data.table(draws@fixef)
  spline1_vars <- names(draws)[grepl('degree = 1[)]1',names(draws))]
  for(v in spline1_vars) draws[, (v) := get(v)/(knot - knot.boundary.left)]
  spline2_vars <- names(draws)[grepl('degree = 1[)]2',names(draws))]
  for(v in spline2_vars) {
    spline1_var <- gsub('degree = 1)2','degree = 1)1',v)
    draws[, (v) := (get(v)-get(spline1_var))/(knot.boundary.right - knot)]
  }
  draws <- draws[, c(spline1_vars,spline2_vars), with=F]
  for(v in c(spline1_vars,spline2_vars)) {
    coefs[variable==v, mean := mean(draws[, get(v)])]
    coefs[variable==v, se := sd(draws[, get(v)])]
  }
  #################################################
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
clean_names <- data.table(variable=c("(Intercept)","as.factor(female)1","as.factor(cohort_group)1940","as.factor(cohort_group)1950",'baseline_cog',"baseline_wealth","baseline_income","mother_edu","child_health",'edu_years',"edu_years_mother","center_age","raceblack",
                                     "edu_years:center_age",
                                     "center_age:raceblack",
                                     "edu_years:raceblack",
                                     "edu_years:center_age:raceblack",
                                     "center_age:edu_years_mother",
                                     "edu_years_mother:raceblack",
                                     "center_age:edu_years_mother:raceblack",
                                     'over_70','over_70:raceblack','center_age:over_70',
                                     "center_age:over_70:raceblack",
                                     'bs(center_age, knots = 4, degree = 1)1',
                                     'bs(center_age, knots = 4, degree = 1)1:raceblack',
                                     'bs(center_age, knots = 4, degree = 1)1:edu_years',
                                     'bs(center_age, knots = 4, degree = 1)2',
                                     'bs(center_age, knots = 4, degree = 1)2:raceblack',
                                     'bs(center_age, knots = 4, degree = 1)2:edu_years',
                                     "mar_rsq","cond_rsq","bic","N","people"),
                          name=c('Intercept','Female','1940 cohort','1950 cohort','Cognition (baseline)','Wealth (baseline)','Income (baseline)','Maternal education','Child health','Years of education','Mother years of education','Age','Black','Age*Edu','Age*Black','Edu*Black','Age*Edu*Black','Age*MEdu','MEdu*Black','Age*MEdu*Black',
                                 'Over70','Over70*Black','Age*Over70','Age*Over70*Black',
                                 'Age (<70)','Age (<70) * Black','Age (<70) * Edu',
                                 'Age (>=70)','Age (>=70) * Black','Age (>=70) * Edu',
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
  hline(j=1:13,i=17,part = 'body', border = fp_border(color='black',width=2)) %>%
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
  print(target = paste0("table2_", file_tag, "_spline.docx"))
