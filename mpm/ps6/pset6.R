# Code for ECO518 MPM Pset 1 


library(tidyverse)
library(readxl)

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps6/")
out <- paste0(dir, "out/")
setwd(dir)

df <- read_xlsx("Guns.xlsx")

