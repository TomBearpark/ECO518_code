# Preliminaries ----
library(tidyverse)
library(sandwich)
library(lmtest)

rm(list = ls(all.names = TRUE))
setwd('/Users/ericqian/Princeton/Courses/ECO 518/homework/hw6')
source('routines.R')

# Settings
nBoot = 500
QQ   = c(.025, .975)



# Question 1 ----
Raw = read.csv('Guns.csv')
df  = Raw
df  = Raw %>% mutate(lvio = log(vio), i_stateid = factor(stateid))
n   = nrow(df)

controls = "incarc_rate + density + avginc + pop + pb1064 + pw1064 + pm1029 + i_stateid"
form     = paste0('lvio ~ shall +', controls)
lm1      = lm(form, data = df)
se      = coeftest(lm1, vcovHC(lm1, type = "HC1"))['shall', 'Std. Error']

betaHat = lm1$coefficients['shall']
summary(lm1)



# Non-parametric bootstrap
Res_np     = replicate(nBoot, np_boot(form, df))
Beta_np    = sapply(Res_np['coeff',], '[[', 'shall')
Beta_np_sd = sd(Beta_np)
CI_np_ef   = quantile(Beta_np, QQ)  # Efron
SE_np      = sapply(Res_np['SE',], '[[', 'shall')
t_np       = (Beta_np - betaHat) / SE_np
CI_np_perc = betaHat - quantile(t_np,QQ)[c(2,1)] * Beta_np_sd
names(CI_np_perc) = names(CI_np_perc)[c(2,1)]



# Parametric bootstrap (residual)
Beta_para    = replicate(nBoot, para_boot(form, df, lm1))
CI_para_ef   = quantile(Beta_para['shall',], QQ)
Beta_para_sd = sd(Beta_para['shall',])


# Clustered bootstrap
Res_clust        = replicate(nBoot, cluster_boot(form, df, 'stateid'))
beta_clust       = sapply(Res_clust['coeff',], '[[', 'shall')
beta_clust_sd    = sd(beta_clust)
CI_clust_ef      = quantile(beta_clust, QQ)
SE_clust         = sapply(Res_clust['SE',], '[[', 'shall')
t_clust          = (beta_clust - betaHat) / SE_clust
CI_clust_perc       = betaHat - quantile(t_clust,QQ)[c(2,1)] * beta_clust_sd
names(CI_np_perc) = names(CI_clust_perc)[c(2,1)]
CI_clust_perc
CI_clust_ef


