# code for MPM pset 7

rm(list = ls())
library(tidyverse) 
library(readxl)
library(sandwich)
library(stargazer)
theme_set(theme_bw())

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps7/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)


