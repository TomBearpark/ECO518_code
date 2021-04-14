rm(list = ls())
if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, quantreg, janitor, xtable)
theme_set(theme_bw())
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/")
setwd(dir);set.seed(1);theme_set(theme_bw())
out <- paste0(dir, "ps9/out/");dir.create(out, showWarnings = FALSE)

######################################################
# Q1 - Probit Marginal Effects 
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

# Take a bootstrap draw for the MLE parameter, given data and initial beta value
draw_params <- function(X, Y, beta0) {
  N <- length(Y)
  sampler <- sample(seq(1, N), N, replace = TRUE)
  Y <- Y[c(sampler)]
  X <- X[c(sampler),]
  MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, 
               method = "BFGS", hessian = TRUE)
  MLE$par
}

# Run MLE estimator for beta - minimise probit negative Log Likelihood 
MLE <- optim(par = beta0, probit.log.lik, X = X, Y = Y, 
             method = "BFGS", hessian = TRUE)

######################################################
# (i) - Marginal effect at the mean for educ

beta <- MLE$par
Xval <- apply(X, 2, mean) %>% as.matrix() # take means of covariates

# Calculate marginal effect for probit model, at Xval for a given variable
probit_marginal_effect <- function(variable, beta, Xval){
  yhat <- t(beta) %*% Xval %>% drop()
  beta[variable,] * dnorm(yhat)
}

# Take a draw of the marginal effect estimate
boot_ME <- function(i, X, Y, beta0, Xval){
  beta <- draw_params(X, Y, beta0)
  tibble(i, coef = probit_marginal_effect("educ", beta, Xval))
}

# Run calculations
me1    <- probit_marginal_effect("educ", beta, Xval) 
draws1 <- map_dfr(1:B, boot_ME, X = X, Y = Y, beta0 = beta0, Xval = Xval)
sd1    <- sd(draws1$coef)

# Plot distribution of the ME draws 
plot_draws <- function(draws, central_val){
  sd <- sd(draws$coef)
  ggplot(draws) + 
    geom_density(aes(x = coef)) + 
    geom_vline(xintercept = central_val, color = "red") + 
    geom_vline(xintercept = c(central_val -1.96 * sd, central_val +1.96 * sd), 
               color = "red", alpha = .6)
}
plot_draws(draws1, me1)

######################################################
# (ii) - Average Partial effect

ape_probit <- function(variable, X, beta){
  N <- dim(X)[1]
  g <- sum(dnorm(X %*% beta)) / N
  beta[variable,] * g
}

boot_APE <- function(i, X, Y, beta0){
  beta <- draw_params(X, Y, beta0)
  tibble(i, coef = ape_probit("educ", X, beta))
}
ape2 <- ape_probit('educ', X, beta) 
draws2 <- map_dfr(1:B, boot_APE, X = X, Y = Y, beta0 = beta0)
sd2 <- sd(draws2$coef)
plot_draws(draws2, ape2)

######################################################
# (iii)

APE_discrete <- function(variable, gap, X, beta){
  N <- dim(X)[1]
  X1 <- X0 <- X
  X1[,variable] <- gap
  X0[,variable] <- 0
  1 / N * sum(pnorm(X1 %*% beta) - pnorm(X0 %*% beta))
}

boot_APE_discrete <- function(i, X, Y, beta0, gap, variable){
  beta <- draw_params(X, Y, beta0)
  tibble(i, coef = APE_discrete(variable, gap, X, beta))
}
ape3 <- APE_discrete('kidslt6', gap = 1,X, beta) 
draws3 <- map_dfr(1:B, boot_APE_discrete, 
                  X = X, Y = Y, beta0 = beta0, gap = 1, variable = "kidslt6")
sd3 <- sd(draws3$coef)
plot_draws(draws3, ape3)

######################################################
# format results
tibble(
  Question = c("(i)", "(ii)", "(iii)"),
  Parameter = c("ME educ at Mean", "APE educ", "APE kidslt6"),
  Value = c(me1, ape2, ape3), 
  SE = c(sd1, sd2, sd3)
  ) %>% 
  xtable(digits = 3) %>% 
  print(include.rownames=FALSE)


######################################################
# Q2
######################################################

df2 <- 
  "ps9/heating.csv" %>% 
  read_csv() %>% 
  mutate(one = 1) %>% 
  janitor::clean_names()


######################################################
# Q3 - Quantile Regressions
######################################################

df3 <- read_csv("ps7/nls.csv")
grid <- seq(0.05, 0.95, 0.01)
ggplot(df3, aes(luwe)) + stat_ecdf(geom = "step")

lm(luwe ~ educ + poly(exper, 2) , data = df3)

quant_coef <- function(tau, df, variable = "educ"){
  quant <- rq(luwe ~ educ + poly(exper, 2), tau = tau, data = df)
  tibble(q = tau, coef = coef(quant)[variable])
}

boot_quant <- function(i, df, grid, variable = "educ"){
  df <- slice_sample(df, prop = 1, replace = TRUE)
  map_dfr(grid, quant_coef, df = df, variable = variable) %>% mutate(draw = i)
}

# Run quantile regressions
quant_df_u  <- map_dfr(1:1000, boot_quant, df = df3, grid = grid)
quant_df    <- map_dfr(grid, quant_coef, df = df3) 

# Calculate summary stats for plotting
quant_df_sd <- quant_df_u %>% group_by(q) %>% summarise(sd = sd(coef)) %>% ungroup()
plot_df     <- left_join(quant_df, quant_df_sd) %>% 
  mutate(min = coef - 1.96 * sd, max = coef + 1.96 * sd)

ggplot(data = plot_df) + 
  geom_errorbar(aes(x= q, ymin= min, ymax = max), alpha  = .5) + 
  geom_point(aes(x = q, y = coef), color = "red") + 
  ggtitle("Coefficient on educ by quantile, SEs from 1000 Bootstrap Draws")+
  ggsave(paste0(out, "3_educ_quantile.png"), height = 4, width = 8)

ggplot(data = plot_df) + 
  geom_ribbon(aes(x= q, ymin= min, ymax = max), alpha  = .3) + 
  geom_line(aes(x = q, y = coef), color = "red") + 
  ggtitle("Coefficient on educ by quantile, SEs from 1000 Bootstrap Draws")
  # ggsave(paste0(out, "3_educ_quantile.png"), height = 4, width = 8)

# Extra plot - a few of the quantiles on a raw scatter plot
df3$q0.95 <- predict(rq(luwe ~ educ + exper, tau = 0.95, data = df3))
df3$q0.05 <- predict(rq(luwe ~ educ + exper, tau = 0.05, data = df3))
ggplot(df3) + 
  geom_point(aes(x = educ, y = luwe)) + 
  geom_line(aes(x = educ, y = q0.05), color = "red") + 
  geom_line(aes(x = educ, y = q0.95), color = "red")
