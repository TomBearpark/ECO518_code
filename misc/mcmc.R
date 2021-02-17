# # MCMC

##########################
# Metropolis
##########################

# target kernel - trivial example. Want to draw from a normal

l <- function(x_i, sigma1 = 1, mu_1 = 0) 
  (1 / (sigma1 * sqrt(2 * pi))) * exp(-0.5 * ((x_i - mu_1) / sigma1)^2)

x <- seq(-5,5,.1)
plot(x, l(x))

metrop <- function(n, eps, l, intial = 0) {
  vec <- vector("numeric", n)
  x <- intial
  vec[0] <- x
  for (i in 2:n) {
    innovation   <- rnorm(1, 0, eps)
    candidate    <- x + innovation
    accept_ratio <- l(candidate) / l(x)
    u <- runif(1, 0, 1)
    if(accept_ratio > u) x <- candidate
    vec[i] <- x
  }
  vec
}

draws <- metrop(10000, 0.5, l)
hist(draws[1:(length(draws) / 2)])
hist(draws[(length(draws) / 2):length(draws)])
hist(draws)
plot(ts(draws))
plot(ts(draws[8000:8100]))
plot(ts(draws[1:300]))
acf(draws)
qqnorm(draws)

##########################
# Gibbs Sampler
##########################
# we want to sample from a multimariate distribution
# situation where we can sample for p(theta1|theta2) and p(theta2|theta1)

# Algorithm:
# Start with an initial theta1 and theta2
# draw for theta1, conditional on theta2 initial value
# then draw for theta2, conditional on the value of theta1 we just had 

# Example - joint normal 
# real population:
library(MASS)
library(dplyr)
library(ggplot2)

# p = 1 - covariance across normals
# theta ~ N(c(0, 0), mat(1,1,1,1))
p <- 0.5
pop <- 
  mvrnorm(n = 100000, mu = c(0,0), Sigma = matrix(c(1, p, 1, p), 2, 2)) %>% 
  as.data.frame()
names(pop) <- c("x", "y")
ggplot(pop) + 
  geom_density_2d_filled(aes(x = x, y = y))

# From results on multivariate normal, we know that 
# theta1|theta2 ~ N(p theta2, 1-p^2)

gibbs <- function(n, initial_x){
  vec <- data.frame(x = initial_x, y = 0) %>% 
    bind_rows(data.frame(x = rep(0, n-1), y = rep(0, n-1)))
  for (i in 2:n){
    vec$x[i] <- rnorm(1, 0.5 * vec$y[i-1], 1 - 0.5^2)
    vec$y[i] <- rnorm(1, 0.5 * vec$x[i], 1 - 0.5^2)
  }
  vec
}
draws <- gibbs(10000, 1, 1)
ggplot(draws) + 
  geom_density_2d_filled(aes(x = x, y = y))





