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
pacman::p_load(IDex2019, tidyverse, ggfortify, patchwork, zoo)
theme_set(theme_bw())
# Load data, and plot a simple time series
load("logcovidcases.RData")
df <- cfd19ts

# Note - this is clearly not stationary. Neither are the differences
autoplot(df) /
autoplot(diff(df)) /
autoplot(diff(df, differences = 2))

######################################
# 1. Fit a 9th order linear AR model to the data. 
######################################

# Estimation of the model can be done with a call to rfvar3(). 
# By default rfvar3() uses an improper prior that shrinks toward 
# persistence. You can omit that prior by setting lambda=NULL, 
# mu=NULL. You can choose whether or not to use the prior.

?rfvar3
AR9 <- rfvar3(
  ydata = df,
  lags = 9,
  xdata = NULL,
  const = TRUE,
  breaks = NULL,
  lambda = NULL,
  mu = NULL,
  ic = NULL,
  sigpar = NULL)

# Check i can get the same thing in base R
baseR_AR <- ar(df, aic = FALSE, order.max = 9, method = "ols")
coefs_df <- data.frame(sims = AR9$By[,,seq(1, dim(AR9$By)[3])], 
                       baseR = baseR_AR$ar)

######################################
# 2. Forecast the next 270 days, using the estimates
######################################
# The y0 agument of the forecasting commands should be dimensioned as 
# a 9 x 1 matrix. It should also be a time series object with start 
# date and frequency. Check that this is true with str(y0).
# If not,use y0 <-ts(matrix(y0, ncol=1), end=c(2021,42), freq=365).
# The c(2021, 42) specifies day 42 (i.e.February11) of 2021.

# Format initial condition
y0 <- df[(length(df) - 8):length(df)]
y0 <- ts(matrix(y0, ncol=1), 
        end=c(2021,42), 
        freq=365)
str(y0)

# A forecast using a single set of parameter values can be generated with 
# fcast()

f <- fcast(y0 = y0,         # values of time series to forecast from
           By = AR9$By,     # fitted AR coefficients from part 1
           Bx = AR9$Bx,     
           xdata = NULL,    # exogenous variables (we don't have any)
           horiz = 270,     # forecast horizon
           const = TRUE, 
           shocks = NULL
           )

# Plot! 
ts.plot(df, xlim = (c(2020.2, 2022)), 
        main = "Historic and forecast USA new Covid-19 Infections", 
        ylab = "New Infections (Log)", xlab = "Date")
points(window(f, start = c(2021,42)), type = "l", col = 2)

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

draws_df %>% 
  pivot_longer(cols = V1:V9, names_to = "lag") %>% 
  mutate(lag = paste0("lag ", substr(lag, 2, 2))) %>% 
  ggplot() + 
  geom_density(aes(x = value)) + 
  facet_wrap(~lag) + 
  ggtitle("Parameter Density Plots")

# A set of draws from the posterior on the forecast can be generated 
# with fcastBand(). This program will make plots showing error bands 
# and will return an array of sampled forecasts, unsorted. 
# It can show uncertainty based on parameter values alone 
# (with whichs=0 or, with whichs left at its default value, 
#   based on both parameter and future shock uncertainty.
?fcastBand
p_vals <- c(10, 32, 50, 68, 90)

f_param_un <- 
  fcastBand(
    pdout = draws,
    y0,
    horiz = 270,
    pctiles = p_vals,
    whichv = NULL,
    whichs = 0,
    main = "Parameter Uncertainty Forecast Bands",
    file = paste0(out, "3_param_only_uncertainty.pdf"), 
    xdata = NULL,
    const = TRUE
  )

fc_draws <- f_param_un$fc %>% 
  as.data.frame.table() %>% 
  mutate(date = rep(1:279, 1000)) %>% 
  group_by(Var3) %>% 
    mutate(draw = cur_group_id()) %>% 
  ungroup() %>% 
  select(-c(Var1, var, Var3), value = Freq)

fc_draws %>% filter(draw %in% seq(1,1000)) %>% 
  mutate(draw = factor(draw)) %>% 
  ggplot() + 
  geom_line(aes( x = date, y = value, color = draw), alpha= 0.3) + 
  theme(legend.position = "none") 

fc_draws %>% 
  group_by(date) %>% 
  summarise(x = quantile(value, p_vals / 100), 
            q = p_vals/ 100) %>% 
  ungroup() %>% 
  mutate(quantile = case_when(
    q %in% c(0.1, 0.9) ~ "10-90", 
    q %in% c(0.32, 0.68) ~ "32-68",
    q == 0.5, "median"
  )) %>% 
  ggplot() + 
  geom_line(aes(x = date, y = x, color = quantile))
p

  



