rm(list = ls())
if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, quantreg)
theme_set(theme_bw())
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/")
setwd(dir);set.seed(1);theme_set(theme_bw())

######################################################
# Q1
######################################################

# Load and format data 
df         <- read_csv("ps8/mroz.csv") %>% mutate(one = 1)
covariates <- c("one","kidslt6", "age", "educ", "nwifeinc")
yvar       <- c("part")
X          <- as.matrix(df[covariates])
Y          <- as.matrix(df[yvar])
N          <- length(Y)
beta0      <- solve(t(X) %*% X, t(X) %*% Y)

B <- 1000 # number of bootstrap draws 

# Negative log-likelihood for MLE
probit.log.lik <- function(beta, X, Y){
  Yhat <- X %*% beta
  p <- pnorm(Yhat)
  -sum((1 - Y) * log(1 - p) + Y * log(p))
}

# Take a bootstrap draw for the MLE parameter
draw_params <- function(X, Y, beta0) {
  N <- length(Y)
  sampler <- sample(seq(1, N), N, replace = TRUE)
  Y <- Y[c(sampler)]
  X <- X[c(sampler),]
  MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, 
               method = "BFGS", hessian = TRUE)
  MLE$par
}

# Compute MLE estimator for beta 
MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, 
             method = "BFGS", hessian = TRUE)

######################################################
# (i)

beta <- MLE$par
Xval <- apply(X, 2, mean) %>% as.matrix() # take means of covariates

# Calculate marginal effect for probit model, at Xval for a given variable
probit_marginal_effect <- function(variable,beta,Xval){
  yhat <- t(beta) %*% Xval %>% drop()
  beta[variable,] * dnorm(yhat)
}

# Take a draw of the marginal effect estimate
boot_ME <- function(i,X,Y, beta0){
  beta <- draw_params(X, Y, beta0)
  tibble(i, coef = probit_marginal_effect("educ", beta, Xval))
}

# Run calculations
me1 <- probit_marginal_effect("educ", beta, Xval) 
draws <- map_dfr(1:B, boot_ME, X = X, Y = Y, beta0 = beta0)
sd1 <- sd(draws$coef)

# Plot distribution of the ME draws 
ggplot(draws) + 
  geom_density(aes(x = coef)) + 
  geom_vline(xintercept = me1, color = "red") + 
  geom_vline(xintercept = c(me1 -1.96 * sd1, me1 +1.96 * sd1), 
             color = "red", alpha = .6)

######################################################
# (ii)
ape_probit <- function(variable, beta)
{
  
}

# (iii)
partial_educ_chunker <- function()
{
  
}
  
# (iv) - bootstrap



######################################################
# Q2
######################################################

######################################################
# Q3
######################################################