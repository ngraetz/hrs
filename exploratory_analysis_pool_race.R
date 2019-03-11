## Read data and reshape to be long by id, age, cohort
library(data.table)
library(lme4)
library(ggplot2)
library(survey)
library(gridExtra)
library(mice)
library(survey)
repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)

hrs <- fread('hrsGform.cl_20190129.csv')
hrs[, cohort := as.Date(as.numeric(birthDate), origin = '1960-01-01')]
hrs[, cohort := as.numeric(substr(cohort,1,4))]
hrs[, id := hhidpn]
hrs <- melt(hrs,
            id.vars = c('id','cohort','femalehrsGform','whitehrsGform','blackhrsGform','hispanichrsGform','otherhrsGform','educYrs'),
            measure.vars = patterns('age_','cogScore_','sampWeight_','wealth_','income_'),
            value.name = c('age','cognitive','pweight','wealth','income'))
hrs[, pweight := as.numeric(pweight)]
hrs[, age := as.numeric(age)]
hrs[, cognitive := as.numeric(cognitive)]
hrs[, wealth := as.numeric(wealth)]
hrs[, income := as.numeric(income)]
hrs[, edu_years := as.numeric(educYrs)]
hrs[whitehrsGform==1, race := 'white']
hrs[blackhrsGform==1, race := 'black']
hrs[hispanichrsGform==1, race := 'hispanic']
hrs[otherhrsGform==1, race := 'other']
hrs <- hrs[age >= 50, ]
hrs[, age_group := cut(age, seq(50,95,5))]
hrs[, age_group := as.numeric(substr(age_group, 2, 3))]
hrs[, cohort_group := cut(cohort, seq(1920,1960,10), dig.lab = 10)]
hrs[, cohort_group := as.numeric(substr(cohort_group, 2, 5))]
hrs[, female := femalehrsGform]
hrs[, grep('hrsGform', names(hrs), value=T) := NULL]
hrs[, variable := NULL]
hrs[female==0, sex := 'male']
hrs[female==1, sex := 'female']
## Drop outliers
for(v in c('wealth','income')) {
  hrs <- hrs[get(v) < quantile(hrs[, get(v)], c(0.99)), ]
  hrs <- hrs[get(v) > quantile(hrs[, get(v)], c(0.01)), ]
}
hrs[, log_income := log(income)]

## Multiple imputation for missing data.
## For some reason, id or pweight messes up MI. It is definitely because it introduces some extreme imbalance (i.e. trying to invert an non-invertible matrix somewhere), but I can't think of why...
## I guess it just over-identifies everythink somewhere so there is zero variation.
sapply(hrs, function(x) sum(is.na(x)))
mi_list <- mice(hrs[, c('cohort','age','cognitive','race','female','edu_years')], m=5)
imp_geo <- mice::complete(mi_list)
imp_geo <- as.data.table(imp_geo)
sapply(imp_geo, function(x) sum(is.na(x)))
imp_hrs <- as.data.table(imp_geo)
hrs[, imp_cognitive := imp_hrs[, cognitive]]
saveRDS(hrs, 'hrs_imputed.RDS')

## Exploratory plots of age trajectories
hrs_design <- svydesign(id=~id, weights=~pweight, data=hrs[!is.na(pweight) & !is.na(cognitive), ])
summarize_cognitive <- function(age, stat) {
  by_race <- function(r,a) {
    ## Subset survey design to this age-race.
    dsub <- subset(hrs_design, race==r & age==a)
    if(stat=='mean') {
      s <- as.data.table(merge(svymean(~imp_cognitive, dsub), confint(svymean(~cognitive, dsub))))
      names(s) <- c('mean','se','lower','upper')
    }
    if(stat=='quantile') {
      s <- as.data.table(svyquantile(~imp_cognitive, dsub, c(.25,.5,.75)))
      names(s) <- c('lower','mean','upper')
    }
    s[, race := r]
    s[, age := a]
    return(s)
  }
  all_race <- rbindlist(lapply(c('white','black'), by_race, a=age))
  return(all_race)
}
all <- rbindlist(lapply(unique(hrs[, age]), summarize_cognitive, stat='mean'))
plot1 <- ggplot() +
  geom_line(data=all,
            aes(y=mean,
                x=age,
                color=race),
            size=2) + 
  geom_ribbon(data=all,
              aes(ymin=lower,
                  ymax=upper,
                  x=age,
                  fill=race),
              alpha=0.5) + 
  ylim(c(0,20)) + 
  labs(x='Age',y='Population mean',title='Survey-weighted mean cognitive score by age (stratified by race)') + 
  theme_minimal()
all <- rbindlist(lapply(unique(hrs[, age]), summarize_cognitive, stat='quantile'))
plot2 <- ggplot() +
  geom_line(data=all,
            aes(y=mean,
                x=age,
                color=race),
            size=2) + 
  geom_ribbon(data=all,
              aes(ymin=lower,
                  ymax=upper,
                  x=age,
                  fill=race),
              alpha=0.5) + 
  ylim(c(0,20)) + 
  labs(x='Age',y='Population distribution',title='Survey-weighted quantiles of cognitive score (0.25-0.75) by age (stratified by race)') + 
  theme_minimal()
grid.arrange(plot1, plot2, nrow=1)

## From Courtney
# 1. Mixed effects model of age trajectories of cognition, stratified by race.
# 2. Mixed effects model of age trajectories of cognition, stratified by baseline cognition (3 groups? <10, 10-25, 25-35).
# Covariates for preliminary models:
# Age (also look at non-linear age terms; likely will need quadratic age term)
# Race
# Gender
# Survey period
# Cohort
# I think we should run models with no SES controls, then including SES controls: education, HH income, HH wealth
library(dplyr)
library(gridExtra)
model_trajectories <- function(group_option, outcome_var, numeric_vars, factor_vars, survey_weight = FALSE, use_REML = FALSE) {
  
  ## Subset data to race or baseline cognition.
  if(grepl('low|middle|high', group_option)) {
    ## Try stratifying by cognition at baseline. Take only those individuals where we observe there first score between ages 50-60.
    model_hrs <- hrs[!is.na(get(outcome_var)), ]
    model_hrs[, age_obs := age]
    model_hrs[, age_obs := min(age_obs, na.rm = T), by = id]
    model_hrs[age==age_obs, baseline_cog := get(outcome_var)]
    model_hrs[, baseline_cog := max(baseline_cog, na.rm = T), by = id] ## Repeat within individual
    model_hrs <- model_hrs[age_obs < 60, ]
    ## Set up baseline groups (<10, 10-25, 25-35).
    model_hrs[baseline_cog<15, baseline := 'low']
    model_hrs[baseline_cog %in% 15:19, baseline := 'middle']
    model_hrs[baseline_cog %in% 20:35, baseline := 'high']
    if(grepl('low', group_option)) model_hrs <- model_hrs[baseline=='low' & !is.na(get(outcome_var)), ]
    if(grepl('middle', group_option)) model_hrs <- model_hrs[baseline=='middle' & !is.na(get(outcome_var)), ]
    if(grepl('high', group_option)) model_hrs <- model_hrs[baseline=='high' & !is.na(get(outcome_var)), ]
    if(grepl('white', group_option)) model_hrs <- model_hrs[race=='white', ]
    if(grepl('black', group_option)) model_hrs <- model_hrs[race=='black', ]
    if(grepl('hispanic', group_option)) model_hrs <- model_hrs[race=='hispanic', ]
  }
  if(group_option %in% c('white','black','hispanic')) model_hrs <- hrs[race==group_option & !is.na(get(outcome_var)), ]
  if(group_option == 'all') model_hrs <- copy(hrs)
  
  ## Subset to respondents with at least 3 responses on outcome.
  id_nums <- model_hrs[, .N , by = id]
  dim(id_nums[N == 1, ])
  model_hrs <- model_hrs[id %in% id_nums[N > 3, id], ]
  model_hrs[, id_factor := as.factor(id)]
  
  ## Rescale all variables to mean/sd=0/1. Keep track of means/sds to transform back after predicting.
  scale_vars <- c(outcome_var, numeric_vars)
  means <- model_hrs[, lapply(.SD, mean, na.rm=TRUE), .SDcols=scale_vars]
  means <- melt(means, measure.vars=scale_vars)
  sds <- model_hrs[, lapply(.SD, sd, na.rm=TRUE), .SDcols=scale_vars]
  sds <- melt(sds, measure.vars=scale_vars)
  setnames(sds, 'value', 'sd')
  setnames(means, 'value', 'mean')
  scales <- merge(means,sds)
  model_hrs <- as.data.table(model_hrs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=scale_vars)) # Subtract mean, divide by SD.
  
  ## If we want age polynomials, we need to calculate orthogonal polynomials. Because our ages are relatively high, age and age^2
  ## are basically perfectly correlated so the model can't converge.
  age_polys <- poly(model_hrs[, age], degree = 2)
  model_hrs[, age_poly_1 := age_polys[,1]]
  model_hrs[, age_poly_2 := age_polys[,2]]
  model_hrs <- as.data.table(model_hrs %>% mutate_at(funs(scale(.) %>% as.vector), .vars=c('age_poly_1','age_poly_2'))) # Subtract mean, divide by SD.
  
  ## Fit model, pull out random effects to calculate individual trajectories. Coef() gives random + fixed effects.
  # model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + poly(age,2) + (poly(age,2)|id_factor)'))
  # model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + (poly(age,2)|id_factor)'))
  # model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + age_poly_1 + age_poly_2 + (age_poly_1 + age_poly_2|id_factor)'))
  model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + age_poly_1*race + age_poly_2*race + (age_poly_1|id_factor)'))
  message(model_formula)
  if(survey_weight) {
    model_hrs <- model_hrs[!is.na(pweight) & pweight > 0, ]
    model_hrs[, pweight := pweight / 10000]
    model1 <- lmer(model_formula, weights=pweight, data=model_hrs, REML=use_REML)
  }
  if(!survey_weight) {
    model1 <- lmer(model_formula, data=model_hrs, REML=use_REML)
  }
  res <- coef(model1)$id_factor
  ids <- rownames(res)
  res <- as.data.table(res)
  res[, id_factor := ids]
  # setnames(res, c('(Intercept)','poly(age, 2)1','poly(age, 2)2'), c('age_int','age_slope1','age_slope2'))
  # setnames(res, c('(Intercept)','age_poly_1','age_poly_2'), c('age_int','age_slope1','age_slope2'))
  # res <- res[, c('id_factor','age_int','age_slope1','age_slope2')]
  setnames(res, c('(Intercept)','racehispanic','raceother','racewhite','age_poly_1','age_poly_2','age_poly_1:racehispanic','age_poly_1:raceother','age_poly_1:racewhite','racehispanic:age_poly_2','raceother:age_poly_2','racewhite:age_poly_2'),
           c('age_int','hispanic_int','other_int','white_int','age_slope1','age_slope2','age_slope1_hispanic','age_slope1_other','age_slope1_white',
             'age_slope2_hispanic','age_slope2_other','age_slope2_white'))
  res <- res[, c('id_factor','age_int','hispanic_int','other_int','white_int','age_slope1','age_slope2','age_slope1_hispanic','age_slope1_other','age_slope1_white','age_slope2_hispanic','age_slope2_other','age_slope2_white')]
  plot_data <- merge(model_hrs, res, by='id_factor')
  # plot_data[, growth_curve := age_int + (age_poly_1 * age_slope1) + (age_poly_2 * age_slope2)]
  plot_data[race=='black', growth_curve := age_int + (age_poly_1 * (age_slope1)) + (age_poly_2 * (age_slope2))]
  plot_data[race=='white', growth_curve := (age_int + white_int) + (age_poly_1 * (age_slope1 + age_slope1_white)) + (age_poly_2 * (age_slope2 + age_slope2_white))]
  plot_data[race=='other', growth_curve := (age_int + other_int) + (age_poly_1 * (age_slope1 + age_slope1_other)) + (age_poly_2 * (age_slope2 + age_slope2_other))]
  plot_data[race=='hispanic', growth_curve := (age_int + hispanic_int) + (age_poly_1 * (age_slope1 + age_slope1_hispanic)) + (age_poly_2 * (age_slope2 + age_slope2_hispanic))]
  plot_data[, growth_curve := growth_curve * scales[variable==outcome_var, sd] + scales[variable==outcome_var, mean]]
  plot_data[, group := group_option]
  if(group_option=='all') plot_data[, group := race]
  # setnames(res, c('(Intercept)','age'), c('age_int','age_slope'))
  # res <- res[, c('id_factor','age_int','age_slope')]
  # plot_data <- merge(model_hrs, res, by='id_factor')
  # plot_data[, growth_curve := age_int + (age * age_slope)]
  # plot_data[, growth_curve := growth_curve * scales[variable=='cognitive', sd] + scales[variable=='cognitive', mean]]
  # plot_data[, age := age * scales[variable=='age', sd] + scales[variable=='age', mean]]
  
  ## Calculate population average trajectory.
  # mean_growth <- data.table(age = seq(min(model_hrs$age), max(model_hrs$age), .01),
  #                           mean_int = mean(res$age_int),
  #                           mean_slope = mean(res$age_slope),
  #                           race = race_option)
  # mean_growth[, growth_curve := mean_int + (age * mean_slope)]
  # mean_growth[, growth_curve := growth_curve * scales[variable=='cognitive', sd] + scales[variable=='cognitive', mean]]
  # mean_growth[, age := age * scales[variable=='age', sd] + scales[variable=='age', mean]]
  plot_data[, N := 1]
  mean_growth <- plot_data[, list(N=sum(N), age_poly_1=mean(age_poly_1), age_poly_2=mean(age_poly_2)), by=c('age','race')]
  # mean_growth <- unique(plot_data[, c('age','age_poly_1','age_poly_2')])
  mean_growth[, age_int := mean(res$age_int)]
  mean_growth[, age_slope1 := mean(res$age_slope1)]
  mean_growth[, age_slope2 := mean(res$age_slope2)]
  for(r in c('white','hispanic','other')) {
    mean_growth[race==r, age_int := age_int + mean(res[, get(paste0(r,'_int'))])]
    mean_growth[race==r, age_slope1 := age_slope1 + mean(res[, get(paste0('age_slope1_',r))])]
    mean_growth[race==r, age_slope2 := age_slope2 + mean(res[, get(paste0('age_slope2_',r))])]
  }
  mean_growth[, growth_curve := age_int + (age_poly_1 * age_slope1) + (age_poly_2 * age_slope2)]
  mean_growth[, growth_curve := growth_curve * scales[variable==outcome_var, sd] + scales[variable==outcome_var, mean]]
  mean_growth[, group := group_option]
  if(group_option=='all') mean_growth[, group := race]
  return(list(model1, plot_data, mean_growth))
  
}


for(outcome_var in c('cognitive')) {
  for(REML in c(FALSE)) {
    for(weight in c(FALSE)) {
      for(ses in c(TRUE, FALSE)) {
        
        if(ses) numeric_vars <- c('edu_years','wealth','log_income') 
        if(!ses) numeric_vars <- NULL
        factor_vars <- c('as.factor(female)','as.factor(cohort_group)')
        
        ## All baseline cognition levels
        all_models <- lapply(c('all'), model_trajectories, outcome_var, numeric_vars, factor_vars, survey_weight = weight, use_REML = REML)
        
        ## Plot underlying age trajectories of the outcome (i.e. predictions based only on random age intercepts/slopes, conditional on the rest of the model).
        plot_data <- rbind(all_models[[1]][[2]])
        mean_growth <- rbind(all_models[[1]][[3]])
        plot_data[, unique_ids := uniqueN(id_factor), by = group]
        mean_growth <- merge(mean_growth, unique(plot_data[,c('group','unique_ids')]), by='group')
        gg1 <- ggplot() + 
          geom_line(data = plot_data,
                    aes(x = age,
                        y = growth_curve,
                        group = id_factor),
                    alpha = 0.1) + 
          geom_line(data = mean_growth,
                    aes(x = age,
                        y= growth_curve),
                    color = 'red',
                    size = 3) + 
          xlim(c(50,95)) + 
          ylim(c(0,27)) + 
          labs(x = 'Age',
               y= 'Cognition score') + 
          theme_minimal() + 
          facet_wrap(~group)
        gg2 <- ggplot() + 
          geom_line(data = mean_growth,
                    aes(x = age,
                        y= growth_curve,
                        color = group),
                    size = 3) +
          geom_label(data = mean_growth[age==60, ],
                     aes(x = age,
                         y = growth_curve,
                         label = unique_ids,
                         group=group)) + 
          xlim(c(50,95)) + 
          ylim(c(0,20)) +
          labs(x = 'Age',
               y= 'Cognition score',
               color = 'Race') + 
          theme_minimal()
        
        ## BY BASELINE COGNITION.
        all_models <- lapply(c('low','middle','high'), model_trajectories, outcome_var, numeric_vars, factor_vars, survey_weight = weight, use_REML = REML)
        
        plot_data <- rbind(all_models[[1]][[2]], all_models[[2]][[2]], all_models[[3]][[2]])
        plot_data <- plot_data[race %in% c('white','black')]
        plot_data[, group := paste0(group,'_',race)]
        plot_data[, group := factor(group, levels=c('low_white','middle_white','high_white','low_black','middle_black','high_black'))]
        mean_growth <- rbind(all_models[[1]][[3]], all_models[[2]][[3]], all_models[[3]][[3]])
        mean_growth <- mean_growth[race %in% c('white','black')]
        mean_growth[, group := paste0(group,'_',race)]
        mean_growth[, group := factor(group, levels=c('low_white','middle_white','high_white','low_black','middle_black','high_black'))]
        plot_data[, unique_ids := uniqueN(id_factor), by = group]
        mean_growth <- merge(mean_growth, unique(plot_data[,c('group','unique_ids')]), by='group')
        gg3 <- ggplot() + 
          geom_line(data = plot_data,
                    aes(x = age,
                        y = growth_curve,
                        group = id_factor),
                    alpha = 0.1) + 
          geom_line(data = mean_growth,
                    aes(x = age,
                        y= growth_curve),
                    color = 'red',
                    size = 3) + 
          xlim(c(50,80)) + 
          ylim(c(0,27)) + 
          labs(x = 'Age',
               y= 'Cognition score') + 
          theme_minimal() + 
          facet_wrap(~group)
        gg4 <- ggplot() + 
          geom_line(data = mean_growth,
                    aes(x = age,
                        y= growth_curve,
                        color = group),
                    size = 3) +
          geom_label(data = mean_growth[age==60,],
                     aes(x = age,
                         y = growth_curve,
                         label = unique_ids,
                         group=group)) + 
          scale_color_manual(values=c('#fd8d3c','#f03b20','#bd0026','#41b6c4','#2c7fb8','#253494')) + 
          xlim(c(50,80)) + 
          ylim(c(0,27)) +
          labs(x = 'Age',
               y= 'Cognition score',
               color='Baseline\ncognition') + 
          theme_minimal()
        
        pdf(paste0('output/', outcome_var, '_ses', ses ,'_weight', weight, '_reml', REML, '.pdf'), height=8, width=12)
        print(gg1)
        print(gg2)
        print(gg3)
        print(gg4)
        dev.off()
        
      }
    }
  }
}
