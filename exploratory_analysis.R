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
            id.vars = c('id','cohort','femalehrsGform','whitehrsGform','blackhrsGform','hispanichrsGform','otherhrsGform'),
            measure.vars = patterns('age_','cogScore_','sampWeight_'))
setnames(hrs, c('value1','value2','value3'), c('age','cognitive','pweight'))
hrs[, pweight := as.numeric(pweight)]
hrs[, age := as.numeric(age)]
hrs[, cognitive := as.numeric(cognitive)]
hrs[whitehrsGform==1, race := 'white']
hrs[blackhrsGform==1, race := 'black']
hrs[hispanichrsGform==1, race := 'hispanic']
hrs[otherhrsGform==1, race := 'other']
hrs[, age := cut(age, seq(60,95,5))]
hrs[, age := as.numeric(substr(age, 2, 3))]
hrs <- hrs[!is.na(age), ]
hrs[, female := femalehrsGform]
hrs[, grep('hrsGform', names(hrs), value=T) := NULL]
hrs[, variable := NULL]

## Multiple imputation for missing data.
## For some reason, id or pweight messes up MI. It is definitely because it introduces some extreme imbalance (i.e. trying to invert an non-invertible matrix somewhere), but I can't think of why...
## I guess it just over-identifies everythink somewhere so there is zero variation.
sapply(hrs, function(x) sum(is.na(x)))
mi_list <- mice(hrs[, c('id','cohort','age','cognitive','race','female')], m=5)
imp_geo <- mice::complete(mi_list)
imp_geo <- as.data.table(imp_geo)
sapply(imp_geo, function(x) sum(is.na(x)))
imp_hrs <- as.data.table(imp_geo)
hrs[, imp_cognitive := imp_hrs[, cognitive]]

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
hrs[, age2 := age^2]
model1 <- lmer(cognitive ~ female + cohort + (1+age|id), data=hrs[race=='white', ], REML = FALSE)

