rm(list=ls())

## Options needed for this run.
g_repo <- 'C:/Users/ngraetz/Documents/repos/g_formula/dev'
hrs_repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
local_cores <- 1
use_image <- TRUE
  
setwd(hrs_repo)
library(data.table)
library(tidyverse)
library(nnet)
library(parallel)
library(survey)
source(paste0(g_repo, "/gfit.R"))

## Get wave indicator
hrs <- fread('hrsGform.cl_20190129.csv')
hrs[, cohort := as.Date(as.numeric(birthDate), origin = '1960-01-01')]
hrs[, cohort := as.numeric(substr(cohort,1,4))]
hrs[, id := hhidpn]
hrs <- melt(hrs,
            id.vars = c('id','cohort','femalehrsGform','whitehrsGform','blackhrsGform','hispanichrsGform','otherhrsGform','educYrs'),
            measure.vars = patterns('age_','cogScore_','sampWeight_','intDate_','wealth_','income_'),
            value.name = c('age','cognitive','pweight','intdate','wealth','income'))
hrs[, wave := as.Date(as.numeric(intdate), origin = '1960-01-01')]
hrs[, wave := as.numeric(substr(wave,1,4))]
hrs[, wave_group := cut(wave, seq(1991,2016,2))]
hrs[, wave_group := as.numeric(substr(wave_group, 2, 5))+1]
hrs[, age := as.numeric(age)]

## Load clean imputed data.
DF <- readRDS(paste0(hrs_repo, '/clean_data/df1d_imputed.RDS'))
setnames(DF, 'hhidpn', 'id')
DF[, cohort := as.Date(as.numeric(birth_date), origin = '1960-01-01')]
DF[, cohort := as.numeric(substr(cohort,1,4))]
DF[, cohort_group := cut(cohort, seq(1920,1960,10), dig.lab = 10)]
DF[, cohort_group := as.numeric(substr(cohort_group, 2, 5))]
#DF <- merge(DF, hrs[, c('id','age','wave_group')], by=c('id','age'))
DF[, married := ifelse(marital_stat=='1.married', 1, 0)]

## Define time-varying (tv) and time-constant (tc) variables and lag everything.
tc_vars <- c('female','cohort_group','educ_years','race','migrant')
tv_vars <- c('imp_cognitive','log_income','wealth')
DF <- as.data.table(DF)
DF[, (paste0('l.', tv_vars)) := shift(.SD), by='id', .SDcols=tv_vars]

## Keep only those ids who were first surveyed in 1998 (War Babies).
## Assign year when each participant was first surveyed, and keep only 1998.
# DF <- DF[, original_wave := min(wave_group), by='id']
# DF <- DF[original_wave == 1998, ]
# DF <- DF[!is.na(cohort_group), ] ## Drop 110 observations who are insanely old or young outliers.

keep_ids <- fread("C:/Users/ngraetz/Downloads/vec_hhidpn_keep.csv")
DF <- DF[id %in% keep_ids$x, ]
# DF[, year := wave_group] ## We need to have a "year" variable to move the simulation forward.
DF <- DF[!is.na(race), ] ## Need to drop 10 observations missing race.
dim(DF)

# Build up the longitudinal model base.
baseForm <- as.formula(paste0('~ ', paste(tc_vars, collapse = ' + ')))

# Formula for longitudinal models
formulas <- list(
  update(baseForm, imp_cognitive ~ . + age * (l.imp_cognitive + l.log_income + l.wealth)),
  update(baseForm, wealth ~ . + age * (l.imp_cognitive + l.log_income + l.wealth)),
  update(baseForm, log_income ~ . + age * (l.imp_cognitive + l.log_income + l.wealth)),
  # update(baseForm, bmi ~ . + age * (l.imp_cognitive + l.log_income + l.wealth + l.bmi + l.stroke)),
  # update(baseForm, stroke ~ . + age * (l.imp_cognitive + l.log_income + l.wealth + l.bmi + l.stroke)),
  update(baseForm, died ~ . + age * (l.imp_cognitive + l.log_income + l.wealth))
)
families <- list(gaussian, gaussian, gaussian, quasibinomial) # families for models
functions <- list(svyglm, svyglm, svyglm, svyglm) # functions for results
if(!(length(formulas) == length(functions) & length(formulas) == length(families))) message('FORMULAS, FAMILIES, FUNCTIONS MISALIGNED')

lags <- list(
  ## Lags
    l.imp_cognitive = function(DF) DF %>% mutate(l.imp_cognitive=lag(imp_cognitive)),
    l.wealth = function(DF) DF %>% mutate(l.wealth=lag(wealth)),
    l.log_income = function(DF) DF %>% mutate(l.log_income=lag(log_income))
    # l.bmi = function(DF) DF %>% mutate(l.bmi=lag(bmi)),
    # l.stroke = function(DF) DF %>% mutate(l.stroke=lag(stroke))
  ## Duration-weighted
    # None for now
)

# Here are the natural deterministic and probabilistic rules.
natural_rules <- list(
  # Deterministic
    # Now handled above in the lags 
  # Probabilistic 
    imp_cognitive = function(DF, models, ...) simPredict(DF, models, 1),
    wealth = function(DF, models, ...) simPredict(DF, models, 2),
    log_income = function(DF, models, ...) simPredict(DF, models, 3),
    # bmi = function(DF, models, ...) simPredict(DF, models, 4),
    # stroke = function(DF, models, ...) simPredict(DF, models, 5),
    died = function(DF, models, ...) simPredict(DF, models, 4)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
intervention_rules_good <- list(
  log_income = function(DF, ...) DF %>% mutate(log_income = 11) %>% select(log_income)
)
intervention_rules_bad <- list(
  log_income = function(DF, ...) DF %>% mutate(log_income = 8) %>% select(log_income)
)

# Lastly, we need a set of rules for each specific effect to be simulated. Below are the rules for simulating the direct effect
# using the natural and intervention courses. These leverage simScenario(), which draws stochastically from the course DF provided. 
direct_effect_rules <- list(
  ## DV of interest
  imp_cognitive = function(DF, models, ...) simPredict(DF, models, 1),
  ## Indirect pathways
  wealth = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 2),
  # bmi = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 4),
  # stroke = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 5),
  ## Direct pathway
  log_income = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 3),
  died = function(DF, models, ...) simPredict(DF, models, 4)
)

boots <- 5 # Number of bootstraps, 100 takes a while
replicationSize <- 1 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
set.seed(80085)
setwd(g_repo)

source("./gfit.R")
DF <- DF[pweight!=0, ]
hrs.design <- svydesign(id=~id, weights=~pweight, data=DF)
bootruns <- lapply(1:boots, function(b) {
  
  message(paste0('Bootstrap ', b))
  
  # Sample individuals with replacement, not rows (this is how we propagate model error right now).
  sampleDF_all <- DF %>%
    select(id) %>%
    unique %>%
    sample_frac(replace=TRUE) %>%
    left_join(DF)
  
  ## For testing a single bootstrap, just use the whole dataset so that I can see where the point estimate would be.
  # sampleDF_all <- DF
  
  # Fit the model
  gfitboot_all <- gfit.init.survey(formulas, families, functions, data=sampleDF_all, survey_design=hrs.design)
  
  # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
  mcDF_all <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_all)) 
  
  ## SURVEY CHANGE: Instead of simulating on a weighted dataset, we have to simulate on an unweighted nationally representative dataset. 
  ## The reason is the simScenario() gets all fucked up because weights in the observed data are individual-specific, but simScenario()
  ## draws from the provided courseDF randomly. We can either fix it this way so that random draws are fine, or don't draw at all and just
  ## use the exact values for each individual from the provided courseDF (I still can't think of a reason why this would be a bad thing to do).
  ## We can use sample_frac() with the survey weight to get an unweighted, nationally representative dataset of the same size as the survey dataset.
  ## HRS CHANGES: This is a bit different for HRS design, I think. I think we want to first subset to the first wave (1998), take a nationally 
  ## representative sample of people 50+, then age them through. It's fine to have lots of different ages as long as we predict/subract out mortality.
  mcDF_all <- mcDF_all %>%
    filter(year==1998) 
  mcDF_all <- mcDF_all %>%
    #filter(age==min(age)) %>%
    select(id,pweight) %>%
    unique %>%
    sample_frac(size=1, weight=pweight, replace=TRUE) %>%
    select(id) %>%
    left_join(mcDF_all)
  
  ## HRS CHANGES: Because data are observed every two years, the lagged time-varying variables in our models are all for 2-year lags. This means
  ## we need to progress the simulation in 2-year steps rather than the default 1-year steps. The alternative is treat "age" as a factor variable
  ## in models to model each transition-probability separately. This is ideal when you have irregular jumps in observation and not that many steps,
  ## but would be annoying with this many ages.
  
  # Run the "natural course" rules
  message('NATURAL COURSE')
  natural_course_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all,
                                              year_step = 2, end_year = 20)
  
  # Intervention: good neighborhood (-2) vs. bad neighborhood (2)
  message('INTERVENTION COURSE 1')
  intervention_good <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules_good,
                                          year_step = 2, end_year = 20)
  message('INTERVENTION COURSE 2')
  intervention_bad <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules_bad,
                                         year_step = 2, end_year = 20)
  
  # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
  # CHANGE: do effects on difference between good and bad interventions, rather than between good intervention and natural course.
  message('DIRECT COURSE')
  direct_effect_DF_all <- progressSimulation(mcDF_all, lags, direct_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good,
                                             year_step = 2, end_year = 20)
  
  # Simulate indirect effects [not tested]
  # message('INDIRECT COURSE 1')
  # indirect_effect_pov <- progressSimulation(mcDF_all, lags, pov_indirect_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good,
  #                                           year_step = 2, end_year = 20)
  # message('INDIRECT COURSE 2')
  # indirect_effect_f_absent <- progressSimulation(mcDF_all, lags, f_absent_indirect_effect_rules, gfitboot_all, natural_DF=intervention_bad, intervention_DF=intervention_good,
  #                                                year_step = 2, end_year = 20)
  
  # Return all courses simulated
  list(natural=natural_course_DF_all,
       intervention_good=intervention_good,
       intervention_bad=intervention_bad,
       direct=direct_effect_DF_all)
       # indirect_effect_pov=indirect_effect_pov,
       # indirect_effect_f_absent=indirect_effect_f_absent)
  
})

## Summarize results.
compile_sims_simple <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
      bs[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    as.data.table()
  df[died==1, sim_death_year := min(age, na.rm=T), by='id']
  df[is.na(sim_death_year), sim_death_year := 99999]
  df[, died := as.numeric(ifelse(year >= sim_death_year, 1, 0)), by='id']
  ## DROP PEOPLE WHO DIED
  df <- df[died==0, ]
  df <- df %>%
    group_by(age, sim) %>%
    dplyr::summarize(
      cognitive = mean(imp_cognitive), 
      wealth = mean(wealth),
      income = mean(log_income)
    ) %>%
    select(-sim) %>%
    summarise_all(
      list(
        mean = mean, 
        lwr = function(x) quantile(x, probs=.025),
        upr = function(x) quantile(x, probs=.975))) %>%
    ungroup %>%
    gather("Metric", "Value", -age) %>%
    mutate(measure=gsub("_[A-z ]*", "", Metric)) %>%
    mutate(statistic=gsub("[A-z ]*_", "", Metric)) %>%
    select(-Metric) %>%
    spread("statistic", "Value") %>%
    mutate(type=name)
  return(df)
}
compile_sims_mort <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
    bs[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    as.data.table()
  df[died==1, sim_death_year := min(age, na.rm=T), by='id']
  df[is.na(sim_death_year), sim_death_year := 99999]
  df[, died := as.numeric(ifelse(year >= sim_death_year, 1, 0)), by='id']
  df <- df %>%
    group_by(age, sim) %>%
    dplyr::summarize(
      died = mean(died),
      cognitive = mean(imp_cognitive)
    ) %>%
    select(-sim) %>%
    summarise_all(
      list(
        mean = mean, 
        lwr = function(x) quantile(x, probs=.025),
        upr = function(x) quantile(x, probs=.975))) %>%
    ungroup %>%
    gather("Metric", "Value", -age) %>%
    mutate(measure=gsub("_[A-z ]*", "", Metric)) %>%
    mutate(statistic=gsub("[A-z ]*_", "", Metric)) %>%
    select(-Metric) %>%
    spread("statistic", "Value") %>%
    mutate(type=name)
  return(df)
}
all_runs <- rbindlist(lapply(c('natural','intervention_good','intervention_bad','direct'), compile_sims_simple, bs = bootruns))
all_runs_mort <- rbindlist(lapply(c('natural','intervention_good','intervention_bad','direct'), compile_sims_mort, bs = bootruns))
all_runs <- rbind(all_runs, all_runs_mort[measure=='died', ])

## Get survey-weighted raw data
make_svy_ci <- function(v) {
  if(v=='~died') hrs.design <- svydesign(id=~id, weights=~pweight, data=DF)
  if(v!='~died') hrs.design <- svydesign(id=~id, weights=~pweight, data=DF[died==0, ])
  cis <- svyby(as.formula(v), ~age, hrs.design, svymean, ci=T, na.rm.all=T)
  cis <- as.data.table(cis)
  setnames(cis, gsub('~','',v), 'mean')
  cis[, upr := mean + 1.96*se]
  cis[, lwr := mean - 1.96*se]
  cis[, measure := gsub('~','',v)]
}
all_cis <- rbindlist(lapply(paste0('~', c(tv_vars,'died')), make_svy_ci))
all_cis[, type := 'Observed']
all_cis[measure=='imp_cognitive', measure := 'cognitive']
all_cis[measure=='log_income', measure := 'income']

## Plot 
cov_plot <- all_runs
cov_obs <- all_cis 
good_col <- '#1f78b4'
bad_col <- '#e31a1c'
nat_col <- '#33a02c'
direct_col <- '#b2df8a'
cols <- c('intervention_good' = '#1f78b4', 'intervention_bad' = '#e31a1c', 'Direct effect' = 'black', 'natural' = '#33a02c', 'indirect_effect_pov' = 'black', 'indirect_effect_marriage' = 'black')
## Figure 1. Fitted natural course from g-formula simulation compared to observed data.
fig2_covs <- ggplot(data=cov_plot[type=='natural',], aes(x=age, y=mean)) +
  geom_line(color=nat_col, size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), fill=nat_col, alpha=.5) +
  geom_point(data=cov_obs, shape=21, fill=nat_col, color='black', size=5) + 
  geom_errorbar(data=cov_obs, aes(ymin=lwr,ymax=upr)) + 
  theme_bw() + 
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated course'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')

## Figure 2. Fitted intervention course with direct effect.
fig3_covs <- ggplot(data=cov_plot %>% filter(type %in% c('intervention_good','intervention_bad')), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  scale_fill_manual(values=cols) + 
  scale_color_manual(values=cols) + 
  theme_bw() +
  theme(strip.text.x = element_text(size = 20),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 15)) +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=FALSE, color=guide_legend(title='Simulated course',override.aes = list(size=10))) + labs(x='Child age',y='Population average (proportion or mean)')
