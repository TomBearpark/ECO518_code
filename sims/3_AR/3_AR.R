############################################################################
# 0 Set up environment, load packages and data
############################################################################
rm(list = ls())
set.seed(1)

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
pacman::p_load(IDex2019, tidyverse, ggfortify, patchwork, zoo, tseries)
theme_set(theme_bw())
# Load data, and plot a simple time series
load("logcovidcases.RData")
df <- cfd19ts

# Note - doesn't look very stationary! Neither do the differences
(autoplot(df) + ggtitle("USA log covid cases, rolling weekly average") ) /
autoplot(diff(df))  + ggtitle("First difference")
ggsave(paste0(out, "/0_diagnostic_plot.png"), height = 5, width = 6)

# Save in a convenient format for later plotting
plot_df <- df %>% as.data.frame() %>% 
  mutate(date = row_number() - length(df))

############################################################################
# 1. Fit a 9th order linear AR model to the data. 
############################################################################

# Estimation of the model can be done with a call to rfvar3(). 
# By default rfvar3() uses an improper prior that shrinks toward 
# persistence. You can omit that prior by setting lambda=NULL, 
# mu=NULL. You can choose whether or not to use the prior.

# ?rfvar3
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

# Check I can get the same thing in base R
baseR_AR <- ar(df, aic = FALSE, order.max = 9, method = "ols")
q1_results <- data.frame(sims = AR9$By[,,seq(1, dim(AR9$By)[3])], 
                       baseR = baseR_AR$ar)

# add a version estimated by OLS to make sure i know whats up
estimate_df <- plot_df %>% as_tibble()
for(i in 1:9){
  estimate_df[paste0("lag", i)] = dplyr::lag(as.vector(estimate_df$x), i)
  if(i!=1) f <- paste0(f, " + lag", i)
  else f <- "lag1"
}
lmcoefs <- lm(data = estimate_df, formula(paste0("x ~ ", f)))
q1_results$lm <- coef(lmcoefs[1])[2:10]

# print coefficients to copy into latex
l <- paste0(round(AR9$Bx,4))
for(n in 1:9) 
  l <- paste0(l, " + ",round(q1_results$sims[n], 4), "y_{t-", n, "}")
l

############################################################################
# 2. Forecast the next 270 days, using the estimates
############################################################################
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
f <- window(f, start = c(2021,43))

png(filename = paste0(out, "2_forecast.png"), height = 1000, width = 1200, res=200)
  ts.plot(df, xlim = (c(2020.2, 2022)), 
          main = "Historic and forecast USA new Covid-19 Infections", 
          ylab = "New Infections (Log)", xlab = "Date")
  points(f, type = "l", col = 2)
dev.off()

############################################################################
# 3. 
# Sample 1000 draws from the posterior distribution of AR coefficients and 
# residual variance to generate a forecast with error bands around the 
# forecast that reflects (only) the uncertainty about the parameters of 
# the model. (Make these 90% and 68% bands.)
############################################################################

# Generation of draws from the posterior for 
# the parameters can then be done with a call to postdraw()

# ?postdraw
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

# Plot the parameter pdfs
draws_df %>% 
  pivot_longer(cols = V1:V9, names_to = "lag") %>% 
  mutate(lag = paste0("Lag ", substr(lag, 2, 2), " coefficient")) %>% 
  ggplot() + 
  geom_density(aes(x = value)) + 
  facet_wrap(~lag) + 
  ggtitle("Parameter Density Plots")

ggsave(paste0(out, "/3_parameter_uncertainty.png"), height = 8, width = 10)

# A set of draws from the posterior on the forecast can be generated 
# with fcastBand(). This program will make plots showing error bands 
# and will return an array of sampled forecasts, unsorted. 
# It can show uncertainty based on parameter values alone 
# (with whichs=0 or, with whichs left at its default value, 
#   based on both parameter and future shock uncertainty.

# ?fcastBand

p_vals <- c(5, 16, 50, 84, 95)

fc_draws <- 
  fcastBand(
    pdout = draws,
    y0,
    horiz = 270,
    pctiles = p_vals,
    # whichv = NULL,
    whichs = 0,
    main = "Parameter Uncertainty Forecast Bands",
    file = paste0(out, "3_param_only_uncertainty.pdf"), 
    xdata = NULL,
    const = TRUE
  )

# Function for converting the results of fcastBand to a dataframe
format_fc <- function(fc_draws){
  fc_draws$fc %>% 
    as.data.frame.table() %>% 
    mutate(date = rep(1:279, 1000)) %>% 
    group_by(Var3) %>% 
      mutate(draw = cur_group_id()) %>% 
    ungroup() %>% 
    select(-c(Var1, var, Var3), value = Freq) %>% 
    filter(date > 9) %>% 
    mutate(date = rep(1:270, 1000))
}

# Function to take quantiles across draws at each date, return a nice plot
plot_quantiles <- function(forecast_bands, p_vals, plot_df, title=NULL,
                           file, out){
  forecast_bands %>% 
    group_by(date) %>% 
    summarise(x = quantile(value, p_vals / 100), 
              q = p_vals/ 100) %>% 
    ungroup() %>% 
    mutate(quantile = case_when(
      q %in% c(0.05, 0.95) ~ "5-95 Forecast Band", 
      q %in% c(0.16, 0.84) ~ "16-84 Forecast Band",
      q == 0.5 ~ "Median Forecast"
    )) %>% 
    bind_rows(plot_df %>% mutate(
              q = 0.5001, quantile = "Historic data")) %>% 
    ggplot() + 
    geom_line(aes(x = date, y = x, group = q, color = quantile)) + 
    ylab("Log COVID 19 New Infections") + xlab("Forecast days ahead")+ 
    theme(legend.title = element_blank())
  ggsave(file = paste0(out, file), height = 4, width = 6)
}

fc_draws <- format_fc(fc_draws)
plot_quantiles(fc_draws, p_vals, plot_df, out = out, 
               file = "3_quantiles_param_only.png")

############################################################################
# (4) Generate error bands that include effects both of uncertainty about the 
# model parameters and uncertainty about future disturbance terms.
############################################################################

fc_draws_full <- 
  fcastBand(
    pdout = draws,
    y0,
    horiz = 270,
    pctiles = p_vals,
    main = "Parameter and Shock Uncertainty Forecast Bands",
    file = paste0(out, "3_param_and_shock_uncertainty.pdf"), 
    xdata = NULL,
    const = TRUE
  )
fc_draws_full <- format_fc(fc_draws_full)
plot_quantiles(fc_draws_full, p_vals, plot_df, out = out, 
               file = "4_quantiles_full_uncertainty.png")

############################################################################
# (5) Plot, on a single graph, 20 of the randomly drawn forecast time 
# series that include the effects of both kinds of uncertainty. 
# (This is a different way of visualizing forecast uncertainty.)
############################################################################
sample_index <- sample(1:1000, 20)

p <- fc_draws_full %>% 
  filter(draw %in% sample_index)  %>% 
  ggplot() + 
  geom_line(aes(x = date, y = value, color = as.factor(draw))) + 
  theme(legend.position = "none") + 
  xlab("Forecast days ahead") + ylab("Log COVID 19 New Infections") 
  
ggsave(p, file = paste0(out, "5_20_random_forecasts.png"), height = 5, width = 5)

q <- fc_draws %>% 
  filter(draw %in% sample_index)  %>% 
  ggplot() + 
  geom_line(aes(x = date, y = value, color = as.factor(draw))) + 
  theme(legend.position = "none") + 
  xlab("Forecast days ahead") + ylab("Log COVID 19 New Infections") 

ggsave(q, file = paste0(out, "5_20_random_forecasts_param_uncertainty.png"), 
       height = 5, width = 5)

############################################################################
# (6) Using the same draws from the posterior and future shocks, 
# evaluate the posterior probability that the rate of new infections 
# is smaller at the end of the forecast period than at the start.
############################################################################

fc_draws_full %>% 
  filter(date %in% c(1, 270)) %>% 
  pivot_wider(names_from = date, names_prefix = "date") %>% 
  mutate(increased = ifelse(date1>date270, 1, 0)) %>% 
  summarise(p = mean(increased))

############################################################################
# (7) Using the unconditional joint pdf of the initial conditions implied 
# by the point-estimate of the parameters, calculate how many standard 
# deviations from the process mean is the initial observation.        
# If your point estimates imply an unstable root, you won’t be able to 
# do this part, so start by calculating the roots and checking whether 
# they are all in the stable region.
############################################################################

coef_vec <- 
  AR9$By[,,] %>% 
  unlist() 

roots <- 
  c(1, -coef_vec) %>% 
  as.vector() %>% 
  polyroot() %>% 
  Mod()

roots %>% sort() %>% round(4)

# we can see they are all outside the unit circle! 

# calculate the variance of the process, using the slides (slide 32)
B          <- matrix(rep(0, 9*9), 9, 9)
B[1,]      <- c(AR9$By)
B[2:9,1:8] <- diag(8)
Sigma      <- matrix(rep(0, 9*9), 9, 9)
Sigma[1,1] <- var(AR9$u)


divide <- function(x) 1/x
# Check we get the same thing as before 
eigen(B)$values %>% divide %>% Mod() %>% sort %>% round(4)

# Sigma matrix for companion system, stacked
Sigma_stacked     <- matrix(rep(0, 81), 81,1)
Sigma_stacked[1,] <- var(AR9$u)

# Solve for Var(y) matrix, stacked
R_y_0_stacked <- solve(diag(81) - kronecker(B,B), Sigma_stacked)
R_y_0_stacked[1]

# Calculate distance
mu <- AR9$Bx/sum(c(1, -AR9$By))
(mu - df[1]) / sqrt(R_y_0_stacked[1])


# Alternative, using slide 33, Solve via doubling algorithm
Omegalast <- Sigma
Wlast     <- B
convergent_path <- c()
for (n in 1:100) {
  Omeganext <- Omegalast + Wlast %*% Omegalast %*% t(Wlast)
  Omegalast <- Omeganext
  Wlast     <- Wlast %*% Wlast
  convergent_path <- c(convergent_path, Omegalast[1,1])
}
plot.ts(convergent_path) # note the quick convergence! 
print(sprintf("Var(y) = %9.3f", Omegalast[1,1]))
# we get the same thing as before

############################################################################
# (8) Use a histogram of the residuals and a normal q-q plot of them to
# assess whether the normality assumption is a reasonable approximation.
############################################################################

residuals <- 
  AR9$u %>% 
  as.data.frame() %>% 
  rename(residuals = 1) 

# Histogram of residuals
residuals %>% 
  ggplot() + 
  geom_density(aes(x = residuals))
ggsave(file = paste0(out, "8_residuals_density.png"), height = 5, width = 5)

# QQplot
png(filename = paste0(out, "8_qqplot.png"), height = 1000, width = 1200, res=200)
qqnorm(residuals$residuals)
dev.off()





