repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
hrs <- readRDS('hrs_imputed_v2.RDS')

## Make exploratory tables
hrs[, waves := ifelse(is.na(cognitive), 0, 1)]
hrs[, waves := sum(waves), by=id]
hrs[, baseline_age := min(age, na.rm=T), by=id]
totals <- unique(hrs[!is.na(edu_cat) & !is.na(race), c('id','waves','race','edu_cat'), with=F])
totals[, N := 1]
grand <- totals[, list(N=sum(N)), by='waves']
grand <- dcast(grand, .~waves, value.var='N')
setnames(grand,'.','Variable')
grand[, Variable := 'N']
race_edu_totals <- totals[, list(N=sum(N)), by=c('waves','race','edu_cat')]
race_edu_totals[, Variable := paste0(race,', ',edu_cat)]
race_edu_totals <- dcast(race_edu_totals, Variable~waves, value.var='N')
age <- hrs[, list(baseline_age=mean(baseline_age)), by='waves']
age <- dcast(age, .~waves, value.var='baseline_age')
setnames(age,'.','Variable')
age[, Variable := 'Baseline age']
table1 <- rbind(race_edu_totals,grand,age)

## Make same table with only those where we observe "baseline" cognition
model_hrs <- hrs[!is.na(cognitive), ]
model_hrs[, age_obs := age]
model_hrs[, age_obs := min(age_obs, na.rm = T), by = id]
model_hrs[age==age_obs, baseline_cog := cognitive]
model_hrs[, baseline_cog := max(baseline_cog, na.rm = T), by = id] ## Repeat within individual
model_hrs <- model_hrs[age_obs < 60, ]
model_hrs[, waves := ifelse(is.na(cognitive), 0, 1)]
model_hrs[, waves := sum(waves), by=id]
model_hrs[, baseline_age := min(age, na.rm=T), by=id]
totals <- unique(model_hrs[!is.na(edu_cat) & !is.na(race), c('id','waves','race','edu_cat'), with=F])
totals[, N := 1]
grand <- totals[, list(N=sum(N)), by='waves']
grand <- dcast(grand, .~waves, value.var='N')
setnames(grand,'.','Variable')
grand[, Variable := 'N']
race_edu_totals <- totals[, list(N=sum(N)), by=c('waves','race','edu_cat')]
race_edu_totals[, Variable := paste0(race,', ',edu_cat)]
race_edu_totals <- dcast(race_edu_totals, Variable~waves, value.var='N')
age <- model_hrs[, list(baseline_age=mean(baseline_age)), by='waves']
age <- dcast(age, .~waves, value.var='baseline_age')
setnames(age,'.','Variable')
age[, Variable := 'Baseline age']
table1 <- rbind(race_edu_totals,grand,age)

model_trajectories <- function(group_option, outcome_var, numeric_vars, factor_vars, survey_weight = FALSE, use_REML = FALSE, use_base_cog=FALSE, edu_race=FALSE) {
  
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
  if(group_option == 'all') model_hrs <- copy(hrs[!is.na(get(outcome_var)), ])
  
  ## Calculate baseline cognition and only keep those where we observe cognition at "baseline" (before age 60)
  model_hrs[, age_obs := age]
  model_hrs[, age_obs := min(age_obs, na.rm = T), by = id]
  model_hrs[age==age_obs, baseline_cog := get(outcome_var)]
  model_hrs[, baseline_cog := max(baseline_cog, na.rm = T), by = id] ## Repeat within individual
  model_hrs <- model_hrs[age_obs < 60, ]
  model_hrs[baseline_cog<15, baseline := 'low']
  model_hrs[baseline_cog %in% 15:19, baseline := 'middle']
  model_hrs[baseline_cog %in% 20:35, baseline := 'high']
  
  if(use_base_cog) numeric_vars <- c(numeric_vars, 'baseline_cog')

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
  if(group_option %in% c('all','low','middle','high')) {
    model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + race*(age_poly_1 + age_poly_2) + (age_poly_1|id_factor)'))
    if(edu_race) model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + edu_years*race + race*(age_poly_1 + age_poly_2) + (age_poly_1|id_factor)'))
  }
  if(group_option %in% c('white','black','hispanic')) model_formula <- as.formula(paste0(outcome_var, ' ~ ', paste(c(factor_vars, numeric_vars), collapse=' + '), ' + age_poly_1 + age_poly_2 + (age_poly_1|id_factor)'))
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
  # setnames(res, c('(Intercept)','racehispanic','raceother','racewhite','age_poly_1','age_poly_2','age_poly_1:racehispanic','age_poly_1:raceother','age_poly_1:racewhite','racehispanic:age_poly_2','raceother:age_poly_2','racewhite:age_poly_2'),
  #          c('age_int','hispanic_int','other_int','white_int','age_slope1','age_slope2','age_slope1_hispanic','age_slope1_other','age_slope1_white',
  #            'age_slope2_hispanic','age_slope2_other','age_slope2_white'))
  # res <- res[, c('id_factor','age_int','hispanic_int','other_int','white_int','age_slope1','age_slope2','age_slope1_hispanic','age_slope1_other','age_slope1_white','age_slope2_hispanic','age_slope2_other','age_slope2_white')]
  # plot_data <- merge(model_hrs, res, by='id_factor')
  # # plot_data[, growth_curve := age_int + (age_poly_1 * age_slope1) + (age_poly_2 * age_slope2)]
  # plot_data[race=='black', growth_curve := age_int + (age_poly_1 * (age_slope1)) + (age_poly_2 * (age_slope2))]
  # plot_data[race=='white', growth_curve := (age_int + white_int) + (age_poly_1 * (age_slope1 + age_slope1_white)) + (age_poly_2 * (age_slope2 + age_slope2_white))]
  # plot_data[race=='other', growth_curve := (age_int + other_int) + (age_poly_1 * (age_slope1 + age_slope1_other)) + (age_poly_2 * (age_slope2 + age_slope2_other))]
  # plot_data[race=='hispanic', growth_curve := (age_int + hispanic_int) + (age_poly_1 * (age_slope1 + age_slope1_hispanic)) + (age_poly_2 * (age_slope2 + age_slope2_hispanic))]
  # plot_data[, growth_curve := growth_curve * scales[variable==outcome_var, sd] + scales[variable==outcome_var, mean]]
  # plot_data[, group := group_option]
  # if(group_option=='all') plot_data[, group := race]
  # # setnames(res, c('(Intercept)','age'), c('age_int','age_slope'))
  # # res <- res[, c('id_factor','age_int','age_slope')]
  # # plot_data <- merge(model_hrs, res, by='id_factor')
  # # plot_data[, growth_curve := age_int + (age * age_slope)]
  # # plot_data[, growth_curve := growth_curve * scales[variable=='cognitive', sd] + scales[variable=='cognitive', mean]]
  # # plot_data[, age := age * scales[variable=='age', sd] + scales[variable=='age', mean]]
  # 
  # ## Calculate population average trajectory.
  # # mean_growth <- data.table(age = seq(min(model_hrs$age), max(model_hrs$age), .01),
  # #                           mean_int = mean(res$age_int),
  # #                           mean_slope = mean(res$age_slope),
  # #                           race = race_option)
  # # mean_growth[, growth_curve := mean_int + (age * mean_slope)]
  # # mean_growth[, growth_curve := growth_curve * scales[variable=='cognitive', sd] + scales[variable=='cognitive', mean]]
  # # mean_growth[, age := age * scales[variable=='age', sd] + scales[variable=='age', mean]]
  # plot_data[, N := 1]
  # mean_growth <- plot_data[, list(N=sum(N), age_poly_1=mean(age_poly_1), age_poly_2=mean(age_poly_2)), by=c('age','race')]
  # # mean_growth <- unique(plot_data[, c('age','age_poly_1','age_poly_2')])
  # mean_growth[, age_int := mean(res$age_int)]
  # mean_growth[, age_slope1 := mean(res$age_slope1)]
  # mean_growth[, age_slope2 := mean(res$age_slope2)]
  # for(r in c('white','hispanic','other')) {
  #   mean_growth[race==r, age_int := age_int + mean(res[, get(paste0(r,'_int'))])]
  #   mean_growth[race==r, age_slope1 := age_slope1 + mean(res[, get(paste0('age_slope1_',r))])]
  #   mean_growth[race==r, age_slope2 := age_slope2 + mean(res[, get(paste0('age_slope2_',r))])]
  # }
  # mean_growth[, growth_curve := age_int + (age_poly_1 * age_slope1) + (age_poly_2 * age_slope2)]
  # mean_growth[, growth_curve := growth_curve * scales[variable==outcome_var, sd] + scales[variable==outcome_var, mean]]
  # mean_growth[, group := group_option]
  # if(group_option=='all') mean_growth[, group := race]
  #return(list(model1, plot_data, mean_growth))
  return(list(model1))
  
}

# for(outcome_var in c('imp_cognitive','cognitive')) {
# for(REML in c(FALSE)) {
# for(weight in c(TRUE, FALSE)) {
# for(ses in c(TRUE, FALSE)) {
outcome_var <- 'cognitive'
REML <- FALSE
weight <- FALSE
ses <- FALSE

if(ses) numeric_vars <- c('edu_years','wealth','log_income') 
if(!ses) numeric_vars <- NULL
factor_vars <- c('as.factor(female)','as.factor(cohort_group)')

pull_models <- function(cat) {
  all_models1 <- model_trajectories(cat, outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=FALSE)
  all_models1 <- all_models1[[1]]
  all_models2 <- model_trajectories(cat, outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=TRUE)
  all_models2 <- all_models2[[1]]
  all_models3 <- model_trajectories(cat, outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=FALSE)
  all_models3 <- all_models3[[1]]
  all_models4 <- model_trajectories(cat, outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=TRUE)
  all_models4 <- all_models4[[1]]
  all_models5 <- model_trajectories(cat, outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=FALSE, edu_race = T)[[1]]
  all_models6 <- model_trajectories(cat, outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML, use_base_cog=TRUE, edu_race = T)[[1]]
  all_models <- list(all_models1, all_models2, all_models3, all_models4, all_models5, all_models6)
  return(all_models)
}
all_models <- pull_models('all')

white1 <- model_trajectories('white', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)
white1 <- white1[[1]]
white2 <- model_trajectories('white', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)
white2 <- white2[[1]]
black1 <- model_trajectories('black', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)
black1 <- black1[[1]]
black2 <- model_trajectories('black', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)
black2 <- black2[[1]]
hispanic1 <- model_trajectories('hispanic', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)
hispanic1 <- hispanic1[[1]]
hispanic2 <- model_trajectories('hispanic', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)
hispanic2 <- hispanic2[[1]]
race_models <- list(white1,black1,hispanic1,white2,black2,hispanic2)

# bs1 <- model_trajectories('white', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs2 <- model_trajectories('white', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs3 <- model_trajectories('black', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs4 <- model_trajectories('black', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs5 <- model_trajectories('hispanic', outcome_var, numeric_vars=NULL, factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs6 <- model_trajectories('hispanic', outcome_var, numeric_vars=c('edu_years','wealth','log_income'), factor_vars, survey_weight = TRUE, use_REML = REML)[[1]]
# bs_models <- list(bs1,bs2,bs3,bs4,bs5,bs6)

tab_model(all_models, adjusted = TRUE, dv.labels = c('1','2','3','4','5','6'), show.ci=FALSE, p.style='asterisk', file = 'all_models.html')
tab_model(race_models, adjusted = TRUE, dv.labels = c('White','White','Black','Black','Hispanic','Hispanic'), show.ci=FALSE, p.style='asterisk', file = 'race_models.html')
# tab_model(bs_models, adjusted = TRUE, dv.labels = c('Low','Low','Middle','Middle','High','High'), show.ci=FALSE, p.style='asterisk', file = 'bs_models.html')
