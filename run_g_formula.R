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
DF <- readRDS('./hrs_imputed.RDS')
DF <- merge(DF, hrs[, c('id','age','wave_group')], by=c('id','age'))

## Define time-varying (tv) and time-constant (tc) variables and lag everything.
tc_vars <- c('female','cohort_group','edu_years','race')
tv_vars <- c('imp_cognitive','log_income','wealth')
DF <- as.data.table(DF)
DF[, (paste0('l.', tv_vars)) := shift(.SD), by='id', .SDcols=tv_vars]

## Keep only those ids who were first surveyed in 1998 (War Babies).
## Assign year when each participant was first surveyed, and keep only 1998.
DF <- DF[, original_wave := min(wave_group), by='id']
DF <- DF[original_wave == 1998, ]
DF <- DF[!is.na(cohort_group), ] ## Drop 110 observations who are insanely old or young outliers.
DF[, year := wave_group] ## We need to have a "year" variable to move the simulation forward.
DF <- DF[!is.na(race), ] ## Need to drop 10 observations missing race.
dim(DF)

# Build up the longitudinal model base.
baseForm <- as.formula(paste0('~ ', paste(tc_vars, collapse = ' + '))) 

# Formula for longitudinal models
formulas <- list(
  update(baseForm, imp_cognitive ~ . + age * (l.imp_cognitive + l.log_income + l.wealth)),
  update(baseForm, wealth ~ . + age * (l.imp_cognitive + l.log_income + l.wealth)),
  update(baseForm, log_income ~ . + age * (l.imp_cognitive + l.log_income + l.wealth))
)
families <- list(gaussian, gaussian, gaussian) # families for models
functions <- list(glm, glm, glm) # functions for results
if(!(length(formulas) == length(functions) & length(formulas) == length(families))) message('FORMULAS, FAMILIES, FUNCTIONS MISALIGNED')

lags <- list(
  ## Lags
    l.imp_cognitive = function(DF) DF %>% mutate(l.imp_cognitive=lag(imp_cognitive)),
    l.wealth = function(DF) DF %>% mutate(l.wealth=lag(wealth)),
    l.log_income = function(DF) DF %>% mutate(l.log_income=lag(log_income))
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
    log_income = function(DF, models, ...) simPredict(DF, models, 3)
)

# To calculate direct/indirect effects, we need an intervention to be enforced within the simulated data under the same natural rules. 
intervention_rules <- list(
  log_income = function(DF, ...) DF %>% mutate(log_income = 7) %>% select(log_income)
)

# Lastly, we need a set of rules for each specific effect to be simulated. Below are the rules for simulating the direct effect
# using the natural and intervention courses. These leverage simScenario(), which draws stochastically from the course DF provided. 
direct_effect_rules <- list(
  ## DV of interest
  imp_cognitive = function(DF, models, ...) simPredict(DF, models, 1),
  ## Indirect pathways
  wealth = function(DF, models, natural_DF, ...) simScenario(DF, models, natural_DF, 2),
  ## Direct pathway
  log_income = function(DF, models, natural_DF, intervention_DF, ...) simScenario(DF, models, intervention_DF, 3)
)

boots <- 5 # Number of bootstraps, 100 takes a while
replicationSize <- 5 # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
set.seed(80085)
setwd(g_repo)

source("./gfit.R")
bootruns <- lapply(1:boots, function(b) {
  
  message(paste0('Bootstrap ', b))
  
  # Sample individuals with replacement not rows
  sampleDF_all <- DF %>%
    select(id) %>%
    unique %>%
    sample_frac(replace=TRUE) %>%
    left_join(DF)
  
  ## For testing a single bootstrap, just use the whole dataset so that I can see where the point estimate would be.
  #sampleDF_all <- DF
  
  # Fit the model
  gfitboot_all <- gfit.init(formulas, families, functions, data=sampleDF_all)
  
  # Replicate this bootstrap an arbitrary amount of times to reduce Monte Carlo error (would need to show this converges)
  mcDF_all <- bind_rows(lapply(1:replicationSize, function(i) sampleDF_all)) 
  
  # Run the "natural course" rules
  natural_course_DF_all <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all)
  
  # Run the "intervention course" rules 
  intervention_DF <- progressSimulation(mcDF_all, lags, natural_rules, gfitboot_all, intervention_rules)
  
  # Simulate direct effect drawing stochastically from either the natural or intervention course according to the direct rules 
  direct_effect_DF <- progressSimulation(mcDF_all, lags, direct_effect_rules, gfitboot_all, natural_DF=natural_course_DF_all, intervention_DF=intervention_DF)
  
  # Simulate indirect effects [none for now]

  # Return all courses simulated
  list(natural=natural_course_DF_all,
       intervention=intervention_DF,
       direct=direct_effect_DF)
  
})

## Summarize results.
compile_sims_simple <- function(name, bs) {
  message(name)
  df <- bind_rows(lapply(1:length(bs), function(i){
    bs[[i]][[name]] %>%
      mutate(sim=i)})) %>%
    group_by(age, sim) %>%
    summarize(
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
all_runs <- rbindlist(lapply(c('natural','intervention','direct'), compile_sims_simple, bs = bootruns))

actualDF <- DF %>%
  group_by(age) %>%
  summarize(
    cognitive = mean(imp_cognitive), 
    wealth = mean(wealth),
    income = mean(log_income)
  ) %>%
  gather("measure", "mean", -age) %>%
  mutate(type = 'observed')

## Plot 
cov_plot <- all_runs
cov_obs <- actualDF 
## Figure 1. Fitted natural course from g-formula simulation compared to observed data.
ggplot(data=cov_plot %>% filter(type %in% c('natural')), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  geom_point(data=cov_obs, aes(fill=type), shape=21, color='black', size=5) + 
  scale_fill_manual(values=c('#7fc97f','#beaed4','red','#ffff99')) + 
  theme_classic() +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated\ncourse'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')

## Figure 2. Fitted intervention course with direct effect.
ggplot(data=cov_plot %>% filter(type %in% c('intervention','direct')), aes(x=age, y=mean)) +
  geom_line(aes(color=type), size=1) +
  geom_ribbon(aes(ymin=lwr, ymax=upr, fill=type), alpha=.5) +
  scale_fill_manual(values=c('#7fc97f','#beaed4','red','#ffff99')) + 
  theme_classic() +
  facet_wrap(~measure, scales = 'free_y') +
  guides(fill=guide_legend(title='Simulated\ncourse'), color=FALSE) + labs(x='Child age',y='Population average (proportion or mean)')
