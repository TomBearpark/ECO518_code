# Question 3

library(tidyverse)
library(sandwich)
library(furrr) # parallel functional programming 
library(xtable)

###########################################################
# Functions 
###########################################################

npm_boot <- function(df, i, formula, beta_hat){
  
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  
  beta <- coef(reg)["X"]
  se   <- sqrt(vcovHC(reg, type="HC1")["X", "X"])
  
  T    <- (beta - beta_hat) / se
  
  tibble(i = i, beta = beta, se = se, T = T)
}

resid_boot <- function(df, i, formula, lm1, beta_hat){
  
  df$Y <- predict(lm1) + sample(lm1$residuals, replace = TRUE)
  reg  <- lm(data = df, as.formula(formula))
  
  beta <- coef(reg)["X"]
  se   <- sqrt(vcovHC(reg, type="HC1")["X", "X"])
  
  T    <- (beta - beta_hat) / se
  
  tibble(i = i, beta = beta, se = se, T = T)
}

MC_draw <- function(i, alpha_0, beta_0, sigma_X_0, B = 500){
  
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
  sd_hat   <- sqrt(vcovHC(fit, type="HC1")["X","X"]) 
  
  t        <- abs(coef(fit)["X"] / sd_hat)
  t_tilde  <- abs(sqrt(N) * beta_hat)
  
  # (a) Asymtotically normal critical value:
  Ca <- abs(qnorm((1-.95)/2))
  a  <- as.numeric((t > Ca))
  
  # (b) nonparametric bootstrap critical value
  draws_npm <- map_dfr(seq(1:B), npm_boot, df = df, 
                       formula = "Y ~ X", beta_hat = beta_hat)
  
  Cb <- quantile(abs(draws_npm$T), .95)
  b  <- as.numeric((t > Cb))
  
  # (c) residual bootstrap critical value.
  draws_resid <- map_dfr(seq(1:B), resid_boot, df = df, lm1 = fit, 
                       formula = "Y ~ X", beta_hat = beta_hat)
  Cc <- quantile(abs(draws_resid$T), .95)
  c  <- as.numeric((t > Cc))
  
  
  # (d)
  Cd <- quantile(abs(sqrt(N) * (draws_npm$beta - beta_hat)), .95)
  d  <- as.numeric((t_tilde > Cd))
  
  # (e)
  Ce <- quantile(abs(sqrt(N) * (draws_resid$beta - beta_hat)), .95)
  e  <- as.numeric((t_tilde > Ce))
  
  tibble(draw = i, t = t, t_tilde = t_tilde, 
         a = a, b = b, c = c, d = d, e = e)
}

###########################################################
# Run code


alpha_0 <- 1
beta_0 <- 0
sigma_X_0 <- function(X) 1

MC_draw(i = 1, alpha_0, beta_0, sigma_X_0, B = 1000)

# Parallel compute: 8 workers 
plan(multisession, workers = 8)

# Part (i)
res_i <- future_map_dfr(seq(1,1000), MC_draw, 
               alpha_0 = alpha_0, beta_0 = beta_0, 
               sigma_X_0 = sigma_X_0, B = 500, 
               .options = future_options(seed = TRUE))

sum_i <- summarise_all(res_i, mean) %>% 
  mutate(part = "i") %>% select(-c(draw, t, t_tilde, part)) %>% t() %>% data.frame()

# Part (ii)
sigma_X_0_ii <- function(X) 0.5 * X

res_ii <- future_map_dfr(seq(1,1000), MC_draw, 
                      alpha_0 = alpha_0, beta_0 = beta_0, 
                      sigma_X_0 = sigma_X_0_ii, B = 500, 
                      .options = future_options(seed = TRUE))

sum_ii <- summarise_all(res_ii, mean) %>% 
  mutate(part = "ii") %>% select(-c(draw, t, t_tilde, part)) %>% t() %>% data.frame()

output <- bind_cols(sum_i, sum_ii) 
names(output)= c("i", "ii")
output %>% 
  xtable(digits = 3)


