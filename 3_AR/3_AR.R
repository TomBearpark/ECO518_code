
# 0 Set up environment, load packages
root <- "/Users/tombearpark/Documents/princeton/1st_year/term2/ECO518_Metrics2/" 
dir <- paste0(root, "sims/exercises/3_AR/")
out <- paste0(dir, "out/")

setwd(dir)
if(!require(IDex2019)) install.packages("IDex2019", repos = NULL, type = "source")
if(!require(pacman)) install.packages("pacman")
pacman::p_load(IDex2019, tidyverse)
