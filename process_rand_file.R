#direct and indirect contributions of life course factors to cognitive impairment in the HRS
#nick graetz, jordan weiss, neal marquez(?)
##notes
# r#var  = measure for respondent at wave #
# wave 1 = 1992
# wave 2 = 1993/1994
# wave 3 = 1995/1996
# wave 4 = 1998
# wave 5 = 2000
# wave 6 = 2002
# wave 7 = 2004
# wave 8 = 2006
# wave 9 = 2008
# wave 10 = 2010
# wave 11 = 2012
# wave 12 = 2014
##load packages and functions
#packages
library(car)
library(corrplot)
library(data.table)
library(dplyr)
library(expss)
library(foreign)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(lme4)
library(mice)
library(summarytools)
library(survey)
options("scipen"=999)
rm(list=ls())
repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
##load data
#hrs rand harmonized version P "randhrs1992_2014v2" downloaded 11 june 2018
hrsGform <- fread("hrsGform.cl_20190129.csv")
##clean data
#identifiers and sample weights
hrsGform.cl <- subset(hrsGform, select=c(hhidpn))
hrsGform.cl$hhid <- hrsGform$hhid
hrsGform.cl$pn <- hrsGform$pn
#cohort
hrsGform.cl$cohort <- hrsGform$hacohort
#sample weights
hrsGform.cl[,paste0("sampWeight_", 1:12)] <- hrsGform[,paste0("r",1:12,"wtresp")]
#interview date (SAS FORMAT: days since 1 January 1960)
hrsGform.cl[,paste0("intDate_", 1:12)] <- hrsGform[,paste0("r",1:12,"iwbeg")]
##bith and death dates
#birth date
hrsGform.cl$birthDate <- hrsGform$rabdate
#death_year
hrsGform.cl$deathYear <- ifelse(!is.na(hrsGform$ranyear), hrsGform$ranyear, hrsGform$radyear) #ignoring NDI and interview discrepancies for now
#death_month
hrsGform.cl$deathMonth <- ifelse(!is.na(hrsGform$ranmonth), hrsGform$ranmonth, hrsGform$radmonth)
#death date
hrsGform.cl$deathDate <- as.numeric(difftime(as.Date(sprintf("%s-%s-%s", hrsGform.cl$deathYear, hrsGform.cl$deathMonth, 1)), as.Date("1960-01-01")))
##demographics
#age
hrsGform.cl[,paste0("age_", 1:12)] <- hrsGform[,paste0("r",1:12,"agey_b")]
#female
hrsGform.cl$femalehrsGform <- ifelse(hrsGform$ragender == "2.female", 1, 0)
#hispanic
hrsGform.cl$hispanichrsGform <- ifelse(hrsGform$rahispan == "1.hispanic", 1, 0)
#white
hrsGform.cl$whitehrsGform <- ifelse((hrsGform$raracem == "1.white/caucasian"), 1, 0)
#black
hrsGform.cl$blackhrsGform <- ifelse((hrsGform$raracem == "2.black/african american"), 1, 0)
#other race
hrsGform.cl$otherhrsGform <- ifelse((hrsGform$raracem == "3.other"), 1, 0)
#migrant
hrsGform.cl$migranthrsGform <- ifelse((hrsGform$rabplace == "11.not us/inc us terr"), 1, 0)
#respondent education, category
hrsGform.cl$educCat <- hrsGform$raeduc
#respondent education, years
hrsGform.cl$educYrs <- hrsGform$raedyrs
#maternal education, years
hrsGform.cl$moEducYrs <- hrsGform$rameduc
#paternal education, years
hrsGform.cl$faEducYrs <- hrsGform$rafeduc
#marital status
hrsGform.cl[,paste0("mstat_", 1:12)] <- hrsGform[,paste0("r",1:12,"mstat")]
##economic
#income
hrsGform.cl[,paste0("income_", 1:12)] <- hrsGform[,paste0("h",1:12,"itot")]
#wealth
hrsGform.cl[,paste0("wealth_", 1:12)] <- hrsGform[,paste0("h",1:12,"atotw")]
#medicaid
hrsGform.cl[,paste0("medicaid_", 1:12)] <- (hrsGform[,paste0("r",1:12,"govmd")])
##chronic conditions
#diabetes
hrsGform.cl[,paste0("diab_", 1:12)] <- (hrsGform[,paste0("r",1:12,"diabe")])
#hypertension
hrsGform.cl[,paste0("hibpr_", 1:12)] <- (hrsGform[,paste0("r",1:12,"hibpe")])
#heart problems
hrsGform.cl[,paste0("heart_", 1:12)] <- (hrsGform[,paste0("r",1:12,"hearte")])
#stroke
hrsGform.cl[,paste0("stroke_", 1:12)] <- (hrsGform[,paste0("r",1:12,"stroke")])
#arthritis
hrsGform.cl[,paste0("arth_", 1:12)] <- (hrsGform[,paste0("r",1:12,"arthre")])
##behavioral health
#bmi
hrsGform.cl[,paste0("bmi_", 1:12)] <- (hrsGform[,paste0("r",1:12,"bmi")])
#current smoke
hrsGform.cl[,paste0("smokCur_", 1:12)] <- (hrsGform[,paste0("r",1:12,"smoken")])
#ever smoke
hrsGform.cl[,paste0("smokEv_", 1:12)] <- (hrsGform[,paste0("r",1:12,"smokev")])
#heavy alcohol consumption (3+ DRINKS ON DRINKING DAYS)
hrsGform.cl[,paste0("alcohol_", 3:12)] <- ifelse(hrsGform[,paste0("r",3:12,"drinkn")] >= 3, 1, 0)
#vigorous physical activity (WAVES 1-6: 3+ TIMES PER WEEK; WAVES 712: > 1 PER WEEK)
hrsGform.cl[,paste0("vigAct_", 1:6)] <- (hrsGform[,paste0("r",1:6,"vigact")])
hrsGform.cl[,paste0("vigAct_", 7:12)] <- ifelse(hrsGform[,paste0("r",7:12,"vgactx")] == "1.every day" | hrsGform[,paste0("r",7:12,"vgactx")] == "2.>1 per week", 1, 0)
#moderate physical activity (WAVES 1-6: N/A; WAVES 712: > 1 PER WEEK)
hrsGform.cl[,paste0("modAct_", 7:12)] <- ifelse(hrsGform[,paste0("r",7:12,"mdactx")] == "1.every day" | hrsGform[,paste0("r",7:12,"mdactx")] == "2.>1 per week", 1, 0)
##cognition
#categorized (1 = NORMAL; 2 = CIND; 3 = DEMENTIA)
hrsGform.cl[,paste0("cogCat_", 4:12)] <- (hrsGform[,paste0("cogfunction",seq(1998,2014,2))])
#continuous (RANGE: 0-27; HIGHER SCORES REFLECT BETTER PERFORMANCE)
hrsGform.cl[,paste0("cogScore_", 4:12)] <- (hrsGform[,paste0("cogtot27_imp",seq(1998,2014,2))])
##multiple imputation on all variables
## save hrsGform.cl as .CSV file
write.csv(hrsGform.cl,"~/hrsGform/output/hrsGform.cl_20190129.csv", row.names=FALSE, na="")


repo <- 'C:/Users/ngraetz/Documents/repos/hrs'
setwd(repo)
hrs <- fread('hrsGform.cl_20190129.csv')
hrs[, cohort := as.Date(as.numeric(birthDate), origin = '1960-01-01')]
hrs[, cohort := as.numeric(substr(cohort,1,4))]
hrs[, id := hhidpn]
hrs <- melt(hrs,
            id.vars = c('id','cohort','femalehrsGform','whitehrsGform','blackhrsGform','hispanichrsGform','otherhrsGform','educYrs','educCat'),
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
## Format categorical education
hrs[educCat %in% c('1.lt high-school'), edu_cat := 'less_highschool']
hrs[educCat %in% c('2.ged','3.high-school graduate'), edu_cat := 'highschool']
hrs[educCat %in% c('4.some college'), edu_cat := 'some_college']
hrs[educCat %in% c('5.college and above'), edu_cat := 'college']

## Multiple imputation for missing data.
## For some reason, id or pweight messes up MI. It is definitely because it introduces some extreme imbalance (i.e. trying to invert an non-invertible matrix somewhere), but I can't think of why...
## I guess it just over-identifies everything somewhere so there is zero variation.
sapply(hrs, function(x) sum(is.na(x)))
# mi_list <- mice(hrs[, c('cohort','age','cognitive','race','female','edu_years')], m=5)
# imp_geo <- mice::complete(mi_list)
# imp_geo <- as.data.table(imp_geo)
# sapply(imp_geo, function(x) sum(is.na(x)))
# imp_hrs <- as.data.table(imp_geo)
# hrs[, imp_cognitive := imp_hrs[, cognitive]]
saveRDS(hrs, 'hrs_imputed_v2.RDS')
