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


# re-write:

###########################################################

df$share_food <- df$foodexp / df$income

Y <- df$share_food
X <- df$log_inc

m_nw <- function(X0, X, Y, h){
  K <- function(Xi) dnorm((Xi - X0) / h)
  w <- sapply(X, K) %>% rbind()
  denom <- rowSums(w)
  w <- w / denom
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = K(X0)/denom)
}

CV <- function(hvals, X, Y, mfunc = m_nw){
  CV_out <- data.frame(h = hvals, CV = 0)
  N <- length(Y)
  i <- 1
  for(h in hvals){
    m <- map_dfr(X, mfunc, X, Y, h = h)
    CV_out$CV[i] <- sum(((Y - m$mhat) / (1 - m$wi_Xi))^2) / N
    i <- 1 + i
  }
  CV_out 
}

findMin <- function(CV_results) {
  minCV <- min(CV_results$CV)
  CV_results$h[CV_results$CV == minCV]
}

CV_nw_results <- CV(hvals = seq(0.1, 1, 0.01) , X = X, Y = Y) 
CV_nw_results %>% plot()
CV_h_nw <- findMin(CV_nw_results)
rangeX <- seq(min(X), max(X), length.out = 100)
mhat_nw <- m_nw(rangeX, X, Y, CV_h_nw)

ggplot() + 
  geom_line(data = mhat_nw, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)


###########################################################
m_ll <- function(X0, X, Y, h){
  
  K <- function(Xi) dnorm((Xi - X0) / h)
  
  z <- matrix(c(rep(1, length(X0)), X0), nrow = 2)
  
  mid <- 0
  for(j in 1:length(X)){
    Zj  <- matrix(c(1, X[j]))
    mid <- mid + K(X[j]) * Zj %*% t(Zj)
  }
  mid <- solve(mid)
  
  weight <- function(Xi) t(z) %*% mid %*% (K(Xi) * matrix(c(1, Xi)))
  
  w <- sapply(X, weight) %>% rbind()
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = weight(X0))
}

CV_ll_results <- CV(hvals = seq(1,2, 0.01), X = X, Y = Y, mfunc = m_ll) 
CV_ll_results %>% plot()
CV_h_ll <- findMin(CV_ll_results)
mhat_ll <- map_dfr(X, m_ll, X, Y, h = CV_h_ll)

ggplot() + 
  geom_line(data = mhat_ll, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)

###########################################################
z <- function(x, p){
  z <- c(1)
  for(j in 1:p) z[j+1] <- x^j
  matrix(z)
}

m_poly <- function(X0, X, Y, h){
  
  p <- h
  mid <- 0
  for(j in 1:length(X)){
    Zj  <- z(X[j], p)
    mid <- mid + Zj %*% t(Zj)
  }
  mid <- solve(mid)
  
  weight <- function(Xi) t(z(X0, p)) %*% mid %*%  z(Xi, p)
  
  w <- sapply(X, weight) %>% rbind()
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = weight(X0))
}

pvals <- seq(1,10)
CV_results_poly <- data.frame(h = pvals, CV = 0)

for(p in pvals){
  print(p)
  try(CV_results_poly$CV[p] <- CV(hvals = p, X, Y, mfunc = m_poly))
}
mhat_poly <- map_dfr(X, m_poly, X, Y, h = 2)

ggplot() + 
  geom_line(data = mhat_poly, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)

