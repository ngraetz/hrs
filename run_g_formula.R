rm(list=ls())

## Options needed for this run.
g_repo <- 'C:/Users/ngraetz/Documents/repos/g_formula/dev'
hrs_repo <- ''
local_cores <- 1
use_image <- TRUE
  
setwd(repo)
library(data.table)
library(tidyverse)
library(nnet)
library(parallel)
source(paste0(g_repo, "/gfit.R"))
  
## Load clean imputed data.
DF <- readRDS('./hrs_imputed.RDS')
  
