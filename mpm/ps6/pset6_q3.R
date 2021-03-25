library(tidyverse)
library(sandwich)
library(furrr)

###########################################################
# Question 3
num <- 1000

alpha_0 <- 1
beta_0 <- 0
sigma_X_0 <- function(X) 1


npm_boot <- function(df, i, formula, beta_hat){
  
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  
  beta <- coef(reg)["X"]
  se   <- sqrt(vcovHC(reg, type="HC1")["X", "X"])
  N    <- nrow(draw)
  
  T    <- sqrt(N) * (beta - beta_hat) / se
  
  tibble(i = i, beta = beta, se = se, T = T)
}

resid_boot <- function(df, i, formula, lm1, beta_hat){
  
  df$Y <- predict(lm1) + sample(lm1$residuals, replace = TRUE)
  reg  <- lm(as.formula(formula), data = df)
  
  beta <- coef(reg)["X"]
  se   <- sqrt(vcovHC(reg, type="HC1")["X", "X"])
  N    <- nrow(df)
  T    <- sqrt(N) * (beta - beta_hat) / se
  
  tibble(i = i, beta = beta, se = se, T = T)
}

MC_draw <- function(i, alpha_0, beta_0, sigma_X_0, B = 500){
  message(i)
  # simulate the data
  N       <- 50
  
  X       <- rnorm(N, mean = 0, sd = 2)
  epsilon <- rt(N, df = 5)
  
  df      <- tibble(
                  X = X, 
                  u = sigma_X_0(X) * epsilon
                ) %>% 
                mutate(Y = alpha_0 + beta_0 * X + u)
  
  # Calculate t-stat
  fit      <- lm(data = df, Y ~ X)
  beta_hat <- coef(fit)["X"]
  sd_hat   <- sqrt(vcovHC(fit, type="HC1")[2,2]) 
  t        <- coef(fit)["X"] / sd_hat
  
  QQ <- c(0.025, 0.975)
  
  # (a) Asymtotically normal critical value:
  Ca <- qnorm(QQ, mean = 0, sd = 1)
  a  <- (1- 1*between(t, Ca[1], Ca[2]))
  
  print("done A")
  print(Ca)
  
  # (b) nonparametric bootstrap critical value
  draws_npm <- map_dfr(seq(1:B), npm_boot, df = df, 
                       formula = "Y ~ X", beta_hat = beta_hat)
  Cb <-  sqrt(N) * (beta_hat - sd_hat * quantile(draws_npm$T, QQ)[c(2,1)] / sqrt(N))
  b <- (1- 1*between(t, Cb[1], Cb[2]))
  print("done B")
  print(Cb)
  
  # (c) residual bootstrap critical value.
  draws_resid <- map_dfr(seq(1:B), resid_boot, df = df, lm1 = fit, 
                       formula = "Y ~ X", beta_hat = beta_hat)
  Cc <- sqrt(N) * (beta_hat - sd_hat * quantile(draws_resid$T, QQ)[c(2,1)] / sqrt(N))
  c <- (1- 1*between(t, Cc[1], Cc[2]))
  print("done C")
  print(Cc)
  
  tibble(draw = i, t = t, a = a, b = b, c = c)
}

MC_draw(i = 1, alpha_0, beta_0, sigma_X_0, B = 10)

# Parallel compute: 8 workers 
plan(multisession, workers = 8)

# Part (i)
res_i <- future_map_dfr(seq(1,1000), MC_draw, 
               alpha_0 = alpha_0, beta_0 = beta_0, 
               sigma_X_0 = sigma_X_0, B = 500)

sum_i <- summarise_all(res_i, mean) %>% 
  mutate(part = "i") %>% select(-c(draw, t)) %>% data.frame()

# Part (ii)
sigma_X_0_ii <- function(X) 0.5 * X
res_ii <- future_map_dfr(seq(1,1000), MC_draw, 
                      alpha_0 = alpha_0, beta_0 = beta_0, 
                      sigma_X_0 = sigma_X_0_ii, B = 500)

sum_ii <- summarise_all(res_ii, mean) %>% 
  mutate(part = "ii") %>% select(-c(draw, t)) %>% data.frame()




