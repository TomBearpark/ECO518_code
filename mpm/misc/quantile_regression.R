check <- function(u, tau){
  ifelse(u<=0, (1 - tau)*(-u), tau * u)
}
plot(check(seq(-1, 1, 0.01), 0.25))
plot(check(seq(-1, 1, 0.01), 0.5))
plot(check(seq(-1, 1, 0.01), 0.75))
quantile_error <- function(tau, Y, q){
  sum(check(Y - q, tau)) / length(Y)
}

Y <- rnorm(10000)
quantile(Y, 0.05)
optim(par = mean(Y), quantile_error, Y = Y, tau = 0.05, 
      method = "Brent", lower = min(Y), upper = max(Y))
