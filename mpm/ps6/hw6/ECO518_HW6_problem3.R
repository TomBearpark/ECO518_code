# Settings ----

library(sandwich)
library(lmtest)
library(tidyverse)
rm(list = ls(all.names = TRUE))
setwd('/Users/ericqian/Princeton/Courses/ECO 518/homework/hw6')
source('routines.R')


# DGP settings
alpha0  = 1
beta0   = 0
N       = 50
df      = 5


# Simulation settings
nMC   = 1000
nBoot = 500

runSim(1, 500)


runSim = function(j, nBoot){
  # Draw data ----
  set.seed(j)
  X   = rnorm(N, mean = 0, sd = sqrt(4))
  eps = rt(N, df) 
  
  
  # Part i ----
  u        = eps
  Y1      = alpha0 + beta0 + u
  data1   = as.data.frame(Y1)
  data1$X = X
  form1   = 'Y1~X'    
  lm1     = lm(form1, data1)
  
  # Get EHW
  tN1 = coeftest(lm1, vcovHC(lm1, type = "HC1"))[,"Std. Error"]

  
  # Get boot
  Beta_np1   = replicate(nBoot, np_boot(form1, data1))
  Beta_para1 = replicate(nBoot, para_boot(form1, data1, lm1))
  
  
  # Part ii ----
  u       = sigma0(X) * eps
  Y2      = alpha0 + beta0 + u
  data2   = as.data.frame(Y2)
  data2$X = X
  form2   = 'Y2~X' 
  lm2    = lm(form2, data2)

  # Get EHW
  tN2 = coeftest(lm2, vcovHC(lm2, type = "HC1"))[,"Std. Error"]
  
  # Get boot
  Beta_np2   = replicate(nBoot, np_boot(form2, data2))
  Beta_para2 = replicate(nBoot, para_boot(form2, data2, lm2))
  
  # Add code that stores (a)-(e)

}


# Functions ---

sigma0 = function(X){
  0.5*X
}