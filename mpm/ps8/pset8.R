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
MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, method = "BFGS", hessian = TRUE)
MLE$par %>% xtable()

# Test its similar to canned R package version 
glm(part ~ one + kidslt6 + age + educ + nwifeinc, 
    family = binomial(link = "probit"), data = df) %>% summary()

# Estimate the hessian
h_hat <- function(beta, X, Y){
  
  Yhat <- X %*% beta
  phi  <- pnorm(Yhat)
  Phi  <- dnorm(Yhat)
  N    <- length(Y) 
  k    <- length(beta)
  
  value <- -phi * (
          Y * (phi + Yhat * Phi) / Phi^2 + 
            (1 - Y) * (phi - Yhat * (1 - Phi)) / (1 - Phi^2))
  
  weight <- function(i) value[i] * X[i, ] %*% t(X[i,])
  h <- matrix(0, nrow = k, ncol = k)
  for(i in 1:N) h <- h + weight(i)
  h / N
}

solve(MLE$hessian)%>% diag() %>% sqrt() 
-h_hat(MLE$par, X, Y) %>% diag() %>% sqrt()


# Score
score <- function(beta, X, Y){
  Yhat <- X %*% beta
  phi <- pnorm(Yhat)
  Phi <- dnorm(Yhat)
  value <- phi / (Phi * (1 - Phi))
  weight <- function(i) value[i] * (Y[i] - Phi[i]) * X[i, ]
  s <- matrix(0, nrow = k)
  for(i in 1:N) s <- s + weight(i)
  s
}
score(MLE$par, X, Y)

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

draws %>% group_by(parameter) %>% summarise(sd = sd(value))

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
