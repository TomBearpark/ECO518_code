######################################
# 0 Set up environment, load packages and data
######################################
rm(list = ls())

# Set paths
root <- paste0("/Users/tombearpark/Documents/princeton/",
               "1st_year/term2/ECO518_Metrics2/")
dir <- paste0(root, "sims/exercises/3_AR/")
out <- paste0(dir, "out/")
setwd(dir)

# Load packages
if(!require(IDex2019)) 
  install.packages("IDex2019", repos = NULL, type = "source")
if(!require(pacman)) install.packages("pacman")
pacman::p_load(IDex2019, tidyverse, ggfortify)

# Load data, and plot a simple time series
load("logcovidcases.RData")
df <- cfd19ts
autoplot(df)

######################################
# 1. Fit a 9th order linear AR model to the data. 
######################################

?rfvar3
AR9 <- rfvar3(df, 
              lags = 9, # number of AR coefs to include
              lambda = NULL, mu = NULL, # specify a flat prior
              ic = NULL                 # specify conditioning on initial observations 
              )

coefs <- AR9$By 
coefs %>% unlist() %>% plot()

# Note, this package gives different results to the standard R version
# why is that?
ar(df, aic = FALSE, order.max = 9)

######################################
# 2. Forecast the next 270 days, using the estimates
######################################

y0 <- df[1:9]
y0 <-ts(matrix(y0, ncol=1), end=c(2021,42), freq=365)

f <- fcast(y0 = y0,    # values of time series to forecast from
           By = coefs, # fitted AR coefficients from part 1
           Bx = 0,     # exogenous variables (we don't have any)
           horiz = 270 # forecast horizon
           )

autoplot(f)

######################################
# 3. 
# Sample 1000 draws from the posterior distribution of AR coefficients and 
# residual variance to generate a forecast with error bands around the 
# forecast that reflects (only) the uncertainty about the parameters of 
# the model. (Make these 90% and 68% bands.)
######################################

# Generation of draws from the posterior for 
# the parameters can then be done with a call to postdraw()

?postdraw
draws <- 
  postdraw(
    AR9,               # fitted AR model 
    n = 1000,          # number of draws 
    nosigprior = TRUE  # specify that we don't have a jeffrey's prior
    )

# Convert to a format I know how to work with
draws_df <- 
  draws$By[,,1:9, 1:1000] %>% 
  t() %>% 
  as_tibble()


# A set of draws from the posterior on the forecast can be generated 
# with fcastBand(). This pro- gram will make plots showing error bands 
# and will return an array of sampled forecasts, unsorted. 
# It can show uncertainty based on parameter values alone 
# (with whichs=0 or, with whichs left at its default value, 
#   based on both parameter and future shock uncertainty.




