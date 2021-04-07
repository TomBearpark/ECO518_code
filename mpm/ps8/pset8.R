rm(list = ls())
library(tidyverse) 
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps8/")
out <- paste0(dir, "out/")
setwd(dir);set.seed(1);theme_set(theme_bw())

#############################################################################
# Question 1
#############################################################################

df <- read_csv("mroz.csv")
df$one <- 1

probit.log.lik <- function(beta, X, Y){
  Xhat <- X %*% beta
  sum(Y * log(pnorm(Xhat)) + (1-Y) * log(1 - pnorm(Xhat)))
}

covariates <- c("one","kidslt6", "age", "educ", "nwifeinc")
X <- as.matrix(df[covariates])
Y <- df$part

# Initialise with OLS
beta0 <- solve(t(X) %*% X) %*% t(X) %*% Y

# Compute MLE estimator for beta 
MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y)

# bootstrap SEs
boot <- function(df, covariates, yvar, i){
  df    <- slice_sample(df, prop = 1, replace = TRUE)
  X     <- as.matrix(df[covariates])
  Y     <- df[yvar]
  beta0 <- solve(t(X) %*% X) %*% t(X) %*% y
  MLE   <- optim(par = beta0, probit.log.lik, X = X, Y = Y)
  tibble(parameter = rownames(MLE$par), value = MLE$par[,], i = 1)
}

draws <- map_dfr(1:1000, boot, df = df, 
                 covariates = covariates, yvar = "part")
draws %>% 
  ggplot() + 
  geom_density(aes(x = value)) + 
  facet_wrap(~parameter, scales = "free")

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

# Calculate 2SLS estimator 
estimate_2SLS <- function(Y, X, R){
  
  n <- length(Y)
  S_rx <- t(R) %*% X / n
  S_rr <- t(R) %*% R / n
  S_ry <- t(R) %*% Y / n

  # calculate beta
  beta <- solve(t(S_rx) %*% solve(S_rr) %*% S_rx)  %*% 
                t(S_rx) %*% solve(S_rr) %*% S_ry
  
  # calculate SE
  u_hat <- 
  
  se <- 
  
  return(beta)
}
estimate_2SLS(Y, X, R) %>% drop()
ivreg(logq ~ logp | mixed + stormy, data = df)
