# Code for ECO518 MPM PSet 1 

###########################################################
# 0 Set up, packages 
###########################################################

library(tidyverse) 
library(readxl)
library(sandwich)
theme_set(theme_bw())

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps6/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(123)

###########################################################
# Question 1
###########################################################

# Contents:
# 1.0 Load in data, run the fixed effects regression
# 1.1 Try bootstrapping in all the different ways 
# 1.2 Plot the bootstrap outputs 
# 1.3 Construct CIs using the cluster bootstrap outputs

###########################################################
# 1.0 Load in data, run the fixed effects regression
B <- 1000

df <- read_xlsx("Guns.xlsx") %>% 
  mutate(log_vio = log(vio), state_fe = factor(stateid))

# linear FE regression..
reg1 <- paste0("log_vio ~ shall + incarc_rate + density + avginc + pop ", 
               "+ pb1064 + pw1064 + pm1029 + state_fe")
lm1 <- lm(data = df, formula = as.formula(reg1))
df$resid <- lm1$residuals
summary(lm1)

# Initial plots... 
plot(lm1)
ggplot(data = df, aes(y = resid, x = shall)) + 
  geom_point(col = 'blue') + 
  geom_abline(slope = 0)

ggplot() + 
  geom_density(data = df, aes(x = resid, color = factor(shall))) 

ggplot() + 
  geom_density(data = df, aes(x = resid, color = factor(shall))) + 
  facet_wrap(~state_fe)

##############################################3
# 1.1 Try bootstrapping in all the different ways 

# Non-parametric bootstrap
npm_boot <- function(df, i, formula){
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}
draws_npm <- map_dfr(seq(1:B), npm_boot, df = df, formula = reg1)
sd_npm <- sd(draws_npm$value)

# Try the residual bootstrap... 
resid_boot <- function(df, i , formula, lm){
  df$log_vio <- predict(lm) + sample(lm$residuals)
  reg        <- lm(data = df, as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}
draws_resid <- map_dfr(seq(1:B), resid_boot, df = df, formula = reg1, lm = lm1)
sd_resid <- sd(draws_resid$value)

# Cluster bootstrap...
cluster_boot <- function(df, i, formula){
  
  draw <- df %>% 
    group_nest(state_fe) %>% 
      slice_sample(prop = 1, replace = TRUE) %>% 
    unnest(c(data))
  
  reg  <- lm(data = draw, formula = as.formula(formula))
  se   <- vcovCL(reg, cluster = ~ state_fe)["shall", "shall"] %>% sqrt()
  
  data.frame(i = i, value = coef(reg)["shall"], se = se)
}

draws_cluster <- map_dfr(seq(1:B), cluster_boot, df = df, formula = reg1)
sd_cluster <- sd(draws_cluster$value)

# Wild bootstrap
wild_boot <- function(df, i, formula, lm){
  
  df$omega   <- sample(c(-1,1), nrow(df) ,replace=T)
  df$log_vio <- predict(lm) + df$omega * lm$residuals
  reg        <- lm(data = df, as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
  
}

draws_wild <- map_dfr(seq(1:B), wild_boot, df = df, formula = reg1, lm = lm1)
sd_wild <- sd(draws_wild$value)

# Wild Cluster bootstrap
wild_cluster_boot <- function(df, i, formula, lm){
  
  nclusters <- length(unique(df$state_fe))
  
  draw <- df %>% 
    group_nest(state_fe) %>% 
      mutate(omega = sample(c(-1,1), nclusters, replace=T)) %>% 
    unnest(c(data))
  
  draw$log_vio <- predict(lm) + draw$omega * lm$residuals
  reg          <- lm(data = draw, as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}

draws_wild_cluster <- map_dfr(seq(1:B), wild_cluster_boot, 
                              df = df, formula = reg1, lm = lm1)
sd_wild_cluster <- sd(draws_wild_cluster$value)

###########################################################
# 1.2 Plot the bootstrap outputs 

# Bind all the data together
draws_npm %>% 
  mutate(Bootstrap = paste0("Non-Parametric, sd = ", round(sd_npm, 4))) %>% 
  bind_rows(
    draws_resid %>% 
      mutate(Bootstrap = paste0("Residual, sd = ", round(sd_resid, 4))) 
  ) %>% 
  bind_rows(
    draws_cluster %>% 
      mutate(Bootstrap = paste0("Clustered, sd = ", round(sd_cluster, 4)))
  ) %>% 
  bind_rows(
    draws_wild %>% 
      mutate(Bootstrap = paste0("Wild, sd = ", round(sd_wild, 4)))
  ) %>% 
  bind_rows(
    draws_wild_cluster %>% 
      mutate(Bootstrap = paste0("Wild Cluster, sd = ", round(sd_wild_cluster, 4)))
  ) %>%
  ggplot() + 
    geom_density(aes(x = value, color = Bootstrap)) + 
    ggtitle(paste0(B, " Bootstrap Draws"))
  
ggplot() + 
  geom_density(data = draws_cluster, aes(x = value)) + 
  ggtitle(paste0(B, " Clustered Bootstrap Draws, Mean is: ", 
                 round(mean(draws_cluster$value), 5))) + 
  geom_vline(xintercept = mean(draws_cluster$value) + 1.96 * sd(draws_cluster$value)) + 
  geom_vline(xintercept = mean(draws_cluster$value) - 1.96 * sd(draws_cluster$value)) 
  
# Compute times plot...
if(require(microbenchmark)){
  times <- microbenchmark::microbenchmark(
    resid_boot(df, 1, reg1, lm1), 
    wild_boot(df, 1, reg1, lm1), 
    npm_boot(df, 1, reg1), 
    cluster_boot(df, 1, reg1), 
    wild_cluster_boot(df, 1, reg1, lm1)
  )
  autoplot(times)
}
###########################################################
# 1.3 Construct CIs using the cluster bootstrap outputs

# 1.3.1 Effron CI
betas <- draws_cluster$value
coverage <- .95
sigma_hat <- sd(betas)
N <- nrow(df)

qLow <- quantile(betas, probs = c((1-coverage)/2)) %>% as.numeric()
qHigh <- quantile(betas, probs = c(1-(1-coverage)/2)) %>% as.numeric()

effron_ci <- c(qLow, qHigh)
# 1.3.1 Percentile-T CI
beta_hat <- coef(lm1)["shall"]

draws_cluster <- draws_cluster %>% 
  mutate(T = sqrt(N) * (value - beta_hat) / se)
Ts <- draws_cluster$T

qLowTilde <- quantile(Ts, probs = c((1-coverage)/2)) %>% as.numeric()
qHighTilde <- quantile(Ts, probs = c(1-(1-coverage)/2)) %>% as.numeric()

percentile_ci <- c(beta_hat - sigma_hat * qHighTilde / sqrt(N), 
                   mu - sigma_hat * qLowTilde / sqrt(N))

# Plot the density function, and the CIs 
ggplot() + 
  geom_density(data = draws_cluster, aes(x = value)) + 
  geom_vline(xintercept = effron_ci, color = "red") + 
  geom_vline(xintercept = percentile_ci, color= "blue") + 
  geom_vline(xintercept = mu, color  = "green")

###########################################################
# Question 3
num <- 1000

alpha_0 <- 1
beta_0 <- 0
sigma_X_0 <- function(X) 1


MC_draw <- function(i, alpha_0, beta_0, sigma_X_0){
  # simulate the data
  X <- rnorm(50, mean = 0, sd = 2)
  epsilon <- rt(50, df = 5)
  df <- tibble(
      X = X, 
      u = sigma_X_0(X) * epsilon) %>% 
    mutate(Y = alpha_0 + beta_0 * X + u)
  
  # Calculate t-stat
  fit <- lm(data = df, Y ~ X)
  t   <- coef(fit)["X"] / sqrt(vcovHC(fit, type="HC0")[2,2]) 
  t   <- abs(t)
  tibble(draw = i, t = t)
}
# Asymtotically normal critical value:


MC_draw(i = 1, alpha_0, beta_0, sigma_X_0)


