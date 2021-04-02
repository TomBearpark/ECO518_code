# code for MPM pset 7

rm(list = ls())
library(tidyverse) 
library(furrr)
theme_set(theme_bw())
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps7/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)
plan(multisession, workers = 8)
###########################################################
# Question 1
###########################################################

# 0 Load data, check it out
df <- read_csv(paste0(dir, "engel.csv")) %>% 
  mutate(log_inc = log(income))

ggplot(df) + 
  geom_point(aes(x = foodexp, y = log_inc))
ggplot(df) + 
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

# Helper function for creating vectors of polynomials 
z <- function(x, p){
  z <- c(1)
  for(j in 1:p) z[j+1] <- x^j
  matrix(z)
}

# Returns weight for observation i at X0
w_i <- function(Xi, X0, X, h, type){
  
  if(type == "kernel") {
    
    w_val <- K(Xi, X0, h) / sum(K(X, X0, h))
    
  }else if (type == "local_linear"){
    
    z   <- matrix(c(1, X0))
    Zi  <- matrix(c(1, Xi))
    mid <- matrix(c(0,0,0,0), nrow = 2)
    
    for(j in 1:length(X)){
      Zj  <- matrix(c(1, X[j]))
      mid <- mid + K(X[j], X0, h) * Zj %*% t(Zj)
    }
    
    w_val <- (t(z) %*% solve(mid) %*% (K(Xi, X0, h) * Zi)) %>% as.numeric()
    
  }else if (type == "poly"){
    p <- h
    Zi <- z(X0, p)
    zx <- z(Xi, p)
    Xmat <- c(rep(1, length(X)))
    for(p in 1:h){
      Xmat <- c(Xmat, X^p)
    }
    Xmat <- matrix(Xmat, nrow = length(X))
    w_val <- t(zx) %*% solve(t(Xmat) %*% Xmat) %*% Zi %>% as.numeric()
      
  }else stop("not implemented")
  
  return(w_val)
}


# Get predicted value at X0
m <- function(X0, X, h, Y, type = "kernel"){
  mhat <- 0
  for(i in 1:length(X)){
    mhat <- mhat + Y[i] * w_i(Xi = X[i], X0 = X0, X = X, h = h, type = type)
  }
  # mhat <- sum(Y * w_i(X = X, X0 = X0, h = h, type = type))
  data.frame(x = X0, m = mhat)
}

# Get cross validation MSE, as a function of h
CV <- function(h, X, Y, type = "kernel"){

  N       <- length(X)
  d       <- data.frame(Y = Y, X = X) 
  ms      <- map_dfr(X, m, X = X, h = h, Y = Y, type = type)
  d$resid <- Y - ms$m
  
  d$w <- 0
  for(i in 1:N){
    d$w[i] <-   w_i(X[i], X0 = X[i], X = X, h = h, type = type)
  }
  (1 / N) * sum((d$resid / (1 - d$w))^2)
}

Y <- df$share_food
X <- df$log_inc

hvals <- seq(0.01, 1, 0.01)

CV_results <- data.frame(h = hvals, 
                         CV = map_dbl(hvals, CV, X = X, Y = Y, type = "kernel"))
CV_results <- filter(CV_results, !is.na(CV))
CV_star <- CV_results$h[CV_results$CV == min(CV_results$CV, na.rm = T)]

# Plot!
ggplot(CV_results) + 
  geom_line(aes(x = h, y= CV)) + 
  ggtitle(paste0("h* = ", CV_star))+ 
  ggsave(paste0(out, "/kernel_CV.png"))

range <- seq(min(X), max(X), length.out = 100)
map_dfr(range, m, X = X, h = CV_star, Y = Y) %>% 
  ggplot() + 
  geom_line(aes(x = x, y = m)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) +
  ylim(c(0,1)) + xlim(5.5, 9)+ 
  ggsave(paste0(out, "/kernel_w_scatter.png"))

###########################################################
# iii) Estimate the regression of share of food expenditures on log income 
#   (i.e. the Engel curve) using Local Linear Regression with 
#   a normal kernel in three steps:

# (a) For Local Linear Regression, calculate and plot the cross-validation 
# criterion as a function of h (the bandwidth) on the interval (1, 2] at
# 0.01 increments (e.g. h = 1.01, 1.02, . . . , 1.99, 2).
# (b) Approximate the optimal bandwidth hC V using a grid search over the 
# interval (1, 2] in 0.01 increments.
# (c) Use the optimal bandwidth to estimate the Local Linear Regression of 
# share of food expenditures on log income at each value of X, 
# where X is the log income data. Plot the estimated regression curve 
# and a scatterplot of the data in the same graph.

CV(h = 1.25, X = X, Y = Y, type = "local_linear")

# This computation is very long!! 
CV_results_ll <- 
  data.frame(h = hvals, 
             CV = future_map_dbl(hvals, CV, X = X, Y = Y, type = "local_linear"))

h_star_ll <- CV_results_ll$h[CV_results_ll$CV == min(CV_results_ll$CV, na.rm = T)]
write_csv(CV_results_ll, paste0(out, "local_linear_CV.csv"))

ggplot(CV_results_ll) + 
  geom_line(aes(x = h, y = CV)) + 
  ggsave(paste0(out, "/local_linear_CV.png"))

map_dfr(range, m, X = X, h = h_star_ll, Y = Y, type = "local_linear") %>% 
  ggplot() + 
  geom_line(aes(x = x, y = m)) + 
  ggtitle(paste0("Local Linear Regression, h = ", h_star_ll)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) +
  ylim(c(0,1)) + xlim(5.5, 9) + 
  ggsave(paste0(out, "/local_linear_w_scatter.png"))

###########################################################
# iv) Estimate the regression of share of food expenditures on log income 
# (i.e. the Engel curve) using Polynomial Series Regression in three steps:
# (a) For Polynomial Series Regression, calculate and plot the cross-validation 
# criterion as a function of p (the order of the polynomial) 
# on the grid {1, 2, . . . , 10}.
# (b) Select the polynomial order that minimizes the cross-validation criterion.
# (c) Use the optimal polynomial order to estimate the Polynomial Series 
# Regression of share of food expenditures on log income at each 
# value of X, where X is the log income data. Plot the estimated regression 
# curve and a scatterplot of the data in the same graph.

pvals <- seq(1,10)
CV_results_poly <- data.frame(h = pvals, CV = 0)

CV(p, X, Y, type = "poly")

for(p in pvals){
  print(p)
  try(CV_results_poly$CV[p] <- CV(p, X, Y, type = "poly"))
}

write_csv(CV_results_poly, paste0(out, "poly_CV.csv"))
ggplot(CV_results_poly) + 
  geom_line(aes(x = h, y = CV)) + 
  ggsave(paste0(out, "/poly_CV.png"))

# All the polys on one plot
get_df_poly <- function(h) {
  map_dfr(range, m, X = X, h = h, Y = Y, type = "poly") %>% mutate(h = h)
}
plot_df_poly <- map_dfr(seq(1,4), get_df_poly)

ggplot(plot_df_poly) + 
  geom_line(aes(x = x, y = m, color = factor(h), group =h)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) 

plot_df_poly %>% filter(h == 2) %>% 
  ggplot() + 
    geom_line(aes(x = x, y = m)) + 
    ggtitle(paste0("Polynomial Series Regression, p = ", 2)) + 
    geom_point(data = df, aes(x = log_inc, y = share_food), alpha = 0.2) +
    ylim(c(0,1)) + xlim(5.5, 9) + 
    ggsave(paste0(out, "poly2_w_scatter.png"))

#########################################################
# Question 2
#########################################################

nw <- function(x, X, Y, h, K = dnorm) {
  
  # Arguments
  # x: evaluation points
  # X: vector (size n) with the predictors
  # Y: vector (size n) with the response variable
  # h: bandwidth
  # K: kernel
  
  # Matrix of size n x length(x) (rbind() is called for ensuring a matrix
  # output if x is a scalar)
  Kx <- rbind(sapply(X, function(Xi) K((x - Xi) / h) / h))
  
  # Weights
  W <- Kx / rowSums(Kx) # Column recycling!
  
  # Means at x ("drop" to drop the matrix attributes)
  drop(W %*% Y)
  
}


df2 <- read_csv("nls.csv") %>% 
  rename(Y = luwe, X = educ, R = exper) %>% 
  mutate(R2 = R^2)

reg1 <- lm(df2, formula = "Y ~ X + R + R2") 
summary(reg1)



# ii)

Y <- df2$Y
R <- df2$R
X <- df2$X

# Find an optimal H using CV
hvals <- seq(0.1, 1, 0.1)
CV_results_2 <- data.frame(h = hvals, 
                        CV = map_dbl(hvals, CV, X = R, Y = X, type = "kernel"))
CV_results_2 <- filter(CV_results_2, !is.na(CV))

ggplot(CV_results_2) + 
  geom_line(aes(x = h, y = CV))

CV_star_2 <- CV_results_2$h[CV_results_2$CV == min(CV_results_2$CV, na.rm = T)]

x_grid <- seq(min(R), max(R), length.out = 100)
plot(R, X)
lines(x_grid, nw(x = x_grid, X = R, Y = X, h = 0.7), col = 2)

plot(R, Y)
lines(x_grid, nw(x = x_grid, X = R, Y = Y, h = 0.7), col = 2)


E_X.R <- map_dfr(x_grid, m, X = X, h = CV_star_2, Y = Y) 

ggplot(E_X.R) + 
  geom_line(aes(x = x, y = m)) + 
  geom_point(data = df2, aes(x = R, y = X))




