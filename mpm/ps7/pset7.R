# code for MPM pset 7

rm(list = ls())
library(tidyverse) 
library(readxl)
library(sandwich)
library(stargazer)
library(patchwork)
theme_set(theme_bw())


dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps7/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)

###########################################################
# Question 1
###########################################################

# 0 Load data, check it out
df <- read_csv(paste0(dir, "engel.csv")) %>% 
  mutate(log_inc = log(income))

ggplot(df) + 
  geom_point(aes(x = foodexp, y = log_inc))
p0 <- ggplot(df) + 
  geom_density(aes(x = log_inc))

###########################################################
# 1.(i) 
# Using a normal kernel and the normal reference rule, estimate 
# the density of log income at each value of X, where X is the 
# log income data. Plot the estimated density of log income.

# Normal kernel function
K <- function(Xi, X0, h){
  z <- (Xi - X0) / h
  (2*pi) ^ (-.5) * exp(- (1 / 2) * z^2 ) 
}

# Density function, returns the density at X0, given data X and bandwidth h
f <- function(X0, X, h){
  N <- length(X)
  fhat <- 0
  for(i in 1:N){
    fhat <- fhat + K(Xi = X[i], X0 = X0, h = h)
  }
  data.frame(X = X0, density = (1 / N*h) * fhat)
}

# Calculate h using normal reference rule 
h_star <- function(X){
  N <- length(X)
  sigma <- sd(X)
  1.059 * sigma  / N ^ (0.2)
}

X <- df$log_inc
hist(X)

h <- h_star(X)
density <- map_dfr(X, f, X = X, h = h)
p1 <- ggplot(density) + 
  geom_line(aes( x = X, y = density))
p1


###########################################################
# 1.ii
# Estimate the regression of share of food expenditures on log income 
# (i.e. the Engel curve) using Kernel Regression with a normal 
# kernel in three steps:

df$share_food <- df$foodexp / df$income
ggplot(df) + 
  geom_point(aes(x = log_inc, y = share_food))

# a)
# For Kernel Regression, calculate and plot the cross-validation criterion 
# as a function of h (the bandwidth) on the interval (0,1] in 
# 0.01 increments (e.g. h = 0.01, 0.02, . . . , 0.99, 1).

hvals <- seq(0.01, 1, 0.01)

Y <- df$share_food
X <- df$log_inc

# Returns weight for observation i at X0
w_i <- function(Xi, X0, X, h){
  K(Xi, X0, h) / sum(K(X, X0, h))
}

# Get predicted value at X0
m <- function(X0, X, h, Y){
  data.frame(x = X0, m = sum(Y * w_i(Xi = X, X0 = X0, X = X, h = h)))
}

# Plot to check it looks reasonable for a random h
map_dfr(X, m, X = X, h = 0.1, Y = Y) %>% 
  ggplot() + 
  geom_line(aes(x = x, y = m)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) + 
  ggtitle("Test, h = 0.1")


CV <- function(h, X, Y){

  N       <- length(X)
  d       <- data.frame(Y = Y, X = X) 
  ms      <- map_dfr(X, m, X = X, h = h, Y = Y)
  d$resid <- Y - ms$m
  
  d$w <- 0
  for(i in 1:N){
    d$w[i] <-   w_i(X[i], X0 = X[i], X = X, h = h)
  }
  (1 / N) * sum((d$resid / (1 - d$w))^2)
}

CV_results <- data.frame(h = hvals, CV = map_dbl(hvals, CV, X = X, Y = Y))
CV_results <- filter(CV_results, !is.na(CV))

CV_star <- CV_results$h[CV_results$CV == min(CV_results$CV, na.rm = T)]

# Plot!
ggplot(CV_results) + 
  geom_line(aes(x = h, y= CV)) + 
  ggtitle(paste0("h* = ", CV_star))

map_dfr(X, m, X = X, h = CV_star, Y = Y) %>% 
  ggplot() + 
  geom_line(aes(x = x, y = m)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) 
