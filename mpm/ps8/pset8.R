rm(list = ls())
library(tidyverse) 
library(xtable)
library(sandwich);library(lmtest)
library(estimatr)

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps8/")
out <- paste0(dir, "out/")
setwd(dir);set.seed(1);theme_set(theme_bw())

#############################################################################
# Question 1
#############################################################################

df <- read_csv("mroz.csv") %>% mutate(one = 1)
covariates <- c("one","kidslt6", "age", "educ", "nwifeinc")
X <- as.matrix(df[covariates])
Y <- df$part

probit.log.lik <- function(beta, X, Y){
  Yhat <- X %*% beta
  p <- pnorm(Yhat)
  # negative log-likelihood
  -sum((1 - Y) * log(1 - p) + Y * log(p))
}

# Initialise with OLS
beta0 <- solve(t(X) %*% X) %*% t(X) %*% Y

# Compute MLE estimator for beta 
MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, 
             method = "BFGS", hessian = TRUE)
MLE$par %>% xtable()

# Test its similar to canned R package version 
glm(part ~ one + kidslt6 + age + educ + nwifeinc, 
    family = binomial(link = "probit"), data = df) %>% summary()

beta <- MLE$par

hessian <- function(X, Y, beta){
  yHat <- X %*% beta
  Phi <- pnorm(yHat)
  phi <- dnorm(yHat)
  k <- length(beta)
  N <- length(Y)
  value <- -phi * (
    Y*(phi + yHat * Phi) / Phi^2 + 
      (1 - Y)* (phi - yHat * ( 1 - Phi)) / (1 - Phi)^2
  )
  weight <- function(i) value[i] * X[i, ] %*% t(X[i, ])
  h <- matrix(0, nrow = k, ncol = k)
  for(i in 1:N) h <- h + weight(i)
  h 
}

Hess <- - hessian(X, Y, beta)
H <- Hess %>% solve() %>% diag() %>% sqrt()

H_df <- tibble(parameter = covariates, sd = H, estimator = "Hessian")

score_estimator <- function(X, Y, beta){
  Yhat <- X %*% beta
  Phi <- pnorm(yHat)
  phi <- dnorm(yHat)
  
  N <- length(Y)           # sample size
  k <- length(beta)         # number of coefficients
  
  value <- Y  * (phi / Phi) - ( 1 - Y) * (phi) / (1 - Phi)
  
  weight <- function(i) (value[i] * X[i, ]) %*% t(value[i] * X[i, ])
  
  s <- matrix(0, ncol = k, nrow = k)
  for(i in 1:N) s <- s + weight(i)
  s 
}
Omega <- score_estimator(X, Y, beta)
O <- Omega %>% solve() %>% diag() %>% sqrt() 
O_df <- tibble(parameter = covariates, sd = O, estimator = "Score")

# Sandwich estimator:
S <- solve(H) %*% Omega %*% solve(H) %>%  diag() %>% sqrt()
S_df <- tibble(parameter = covariates, sd = S, estimator = "Sandwich")

# information estimator
I_estimator <- function(X, Y, beta){
  Yhat <- X %*% beta
  Phi <- pnorm(yHat)
  phi <- dnorm(yHat)
  
  N <- length(Y)           # sample size
  k <- length(beta)         # number of coefficients
  
  value <- phi^2 / (Phi * ( 1 - Phi))
  
  weight <- function(i) value[i] * X[i, ] %*% t(X[i, ])
  
  I <- matrix(0, ncol = k, nrow = k)
  for(i in 1:N) I <- I + weight(i)
  I
}

I <- I_estimator(X, Y, beta) %>% solve() %>% diag() %>% sqrt()
I_df <- tibble(parameter = covariates, sd = I, estimator = "I")

# bootstrap SEs
boot <- function(df, covariates, yvar, i){
  df    <- slice_sample(df, prop = 1, replace = TRUE)
  X     <- as.matrix(df[covariates])
  Y     <- df[yvar] %>% as.matrix()
  beta0 <- solve(t(X) %*% X) %*% t(X) %*% Y
  MLE   <- optim(par = beta0, probit.log.lik, X = X, Y = Y)
  tibble(parameter = rownames(MLE$par), value = MLE$par[,], i = 1)
}

draws <- map_dfr(1:1000, boot, df = df, 
                 covariates = covariates, yvar = "part")
draws %>% 
  ggplot() + 
  geom_density(aes(x = value)) + 
  facet_wrap(~parameter, scales = "free")

boot_df <- draws %>% group_by(parameter) %>% summarise(sd = sd(value)) %>% 
  mutate(estimator = "boot")

results <- bind_rows(boot_df, I_df) %>% 
  bind_rows(H_df) %>% bind_rows(O_df) %>% bind_rows(S_df)
pivot_wider(results, id_cols = "estimator", 
            names_from = "parameter", values_from = "sd") %>% xtable(digits = 4)

# iii)
# Find a confidence elipsoid around kidslt6 and educ

#############################################################################
# Question 3
#############################################################################

# Load and format the data
df <- read_csv("fish.csv") %>% 
  mutate(c = 1)
Y <- as.matrix(df[c("logq")])
X <- as.matrix(df[c("c", "logp")])
R <- as.matrix(df[c("c", "mixed","stormy")])

ggplot(data = df) + 
  geom_point(aes(x = logp, y = logq, color = mixed))
ggplot(data = df) + 
  geom_point(aes(x = logp, y = logq, color = stormy))
lm1 <- lm(logp ~ mixed + stormy , data = df) 
first_stage <- coeftest(lm1, vcov = vcovHC(lm1, "HC1")) 

# Calculate 2SLS estimator 
estimate_2SLS <- function(Y, X, R){
  
  N <- length(Y)
  S_rx <- t(R) %*% X / n
  S_rr <- t(R) %*% R / n
  S_ry <- t(R) %*% Y / n
  
  # calculate beta
  beta <- solve(t(S_rx) %*% solve(S_rr) %*% S_rx)  %*% 
                t(S_rx) %*% solve(S_rr) %*% S_ry
  
  # calculate SE
  u_hat <- Y - X %*% beta
  Omega <- diag(N)
  diag(Omega) <- u_hat^2
  
  G <- -1/N * t(R) %*% X
  W <-  1/N * solve(t(R) %*% R)
  O <-  1/N * t(R) %*% Omega %*% R
  
  gmm_var <-  solve(t(G) %*% W %*% G) %*% 
                t(G) %*% W %*% O %*% W %*% G %*% 
              solve(t(G) %*% W %*% G)
  
  sd_gmm <- 1/N * gmm_var %>% diag %>% sqrt()
  
  list(results = tibble(coefficient = drop(beta), sd = sd_gmm), 
       omega = O)
}
results_2sls <- estimate_2SLS(Y, X, R)
results_2sls$results %>% xtable()
# tidy(iv_robust(Y ~ X - 1 | R, data = df)) # compare to canned function

# Use the estimate obtained by 2sls to do two setp GMM

initial <- results_2sls
estimate_2STEP <- function(inital, Y, X, R){
  
  omega <- initial$omega
  W <- solve(omega)
  
  N <- length(Y)
  S_rx <- t(R) %*% X / n
  S_rr <- t(R) %*% R / n
  S_ry <- t(R) %*% Y / n
  
  beta <- solve(t(S_rx) %*% W %*% S_rx) %*% t(S_rx) %*% W %*% S_ry
  
  u_hat <- Y - X %*% beta
  Omega <- diag(N)
  diag(Omega) <- u_hat^2
  
  G <- -1/N * t(R) %*% X
  W <-  1/N * solve(t(R) %*% R)
  O <-  1/N * t(R) %*% Omega %*% R
  
  gmm_var <-  solve(t(G) %*% W %*% G) %*% 
    t(G) %*% W %*% O %*% W %*% G %*% 
    solve(t(G) %*% W %*% G)
  
  sd_gmm <- 1/N * gmm_var %>% diag %>% sqrt()
  
  list(results = tibble(coefficient = drop(beta), sd = sd_gmm), 
       omega = O)
}
results_2STEP <- estimate_2STEP(initial, Y,X,R)
results_2STEP$results %>% xtable
