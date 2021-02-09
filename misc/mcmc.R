# target kernel - trivial example. Want to draw from a normal
l <- function(x_i, sigma1 = 1, mu_1 = 0) 
  (1 / (sigma1 * sqrt(2 * pi))) * exp(-0.5 * ((x_i - mu_1) / sigma1)^2)

x <- seq(-5,5,.1)
plot(x, l(x))

metrop <- function(n, eps, l) {
  vec <- vector("numeric", n)
  x <- 0
  vec[0] <- x
  for (i in 2:n) {
    innovation <- runif(1, -eps, eps)
    candidate <- x + innovation
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
plot(ts(draws))

data.frame(draws = draws) %>% 
  ggplot() + 
  geom_density(aes(x = draws))

plot(ts(draws[8000:8100]))
acf(draws)
qqnorm(draws)
