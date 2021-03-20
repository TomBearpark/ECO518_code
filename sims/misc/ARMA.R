# slide 9 - two MA processes
pacman::p_load(tidyr, ggplot2, dplyr)
df <- 
  tibble(epsilon = rnorm(1000, mean = 0, sd = 2)) %>% 
  mutate(nebla = 0.5 * epsilon, t = row_number())

var(df$epsilon) == var(df$nebla) * 4

df <- df %>% 
  mutate(Ya = epsilon + 0.5 * lag(epsilon), 
         Yb = nebla + 2 * lag(nebla))

df %>% 
  pivot_longer(cols = c(Ya, Yb)) %>% 
  ggplot() + 
  geom_line(aes(x = t, y = value, color = name))

acf(df$Ya[!is.na(df$Ya)])
acf(df$Yb[!is.na(df$Yb)])


# Slide 18, estimating an MA
Y  <- matrix(c(1, 0, -1))
b  <- 0.5

L <- function(b, sigma2, Y){
  message("This function assumes an MA1 process!")
  N <- length(Y)
  mu <- rep(0, length(Y))
  sigma <- diag((1+b^2) * sigma2, nrow = N, ncol = N)
  for (i in 1:N-1){
    sigma[i, i + 1] <- sigma[i+1, i] <- b * sigma2
  }
  # Calculate the likelihood
  -0.5 * log(det(sigma)) - 0.5 * t(Y - mu) %*% solve(sigma) %*% (Y - mu)
}

L(0.5, 1, Y )
# test it -find the MLE
df <- 
  tibble(epsilon = rnorm(100)) %>% 
  mutate(l_epsilon = dplyr::lag(epsilon), 
         y = epsilon + b * l_epsilon) %>% 
  filter(!is.na(y))

Y <- as.matrix(df["y"])
L1 <- function(param) L(b = param[1], sigma2 = 1, Y = Y)
MLE <- optim(c(0.1), L1)











