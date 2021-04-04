rm(list = ls())
library(tidyverse) 

theme_set(theme_bw())
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps8/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)

df <- read_csv("mroz.csv")
df$one <- 1

probit.log.lik <- function(beta, X, Y){
  mu.hat <- X %*% beta
  as.numeric(mu.hat)
  sum(Y * log(pnorm(mu.hat)) + (1-Y) * log(1 - pnorm(mu.hat)))
}

covariates <- c("one","kidslt6", "age", "educ", "nwifeinc")
X <- as.matrix(df[covariates])
Y <- df$part
y <- as.matrix(Y)
beta0 <- solve(t(X) %*% X) %*% t(X) %*% y

optim(par = beta0, probit.log.lik, X = X, Y = Y)

