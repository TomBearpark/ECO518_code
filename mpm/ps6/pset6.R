# Code for ECO518 MPM PSet 1 

###########################################################
# 0 Set up, packages 
###########################################################
rm(list = ls())
library(tidyverse) 
library(readxl)
library(sandwich)
library(stargazer)
theme_set(theme_bw())

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps6/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)

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
N <- nrow(df)

# linear FE regression..
reg1 <- paste0("log_vio ~ shall + incarc_rate + density + avginc + pop ", 
               "+ pb1064 + pw1064 + pm1029 + state_fe")
lm1 <- lm(data = df, formula = as.formula(reg1))
df$resid <- lm1$residuals

ses <- sqrt(diag(vcovCL(lm1, cluster = ~ state_fe)))
stargazer(lm1, keep = "shall", se = list(ses))
beta_hat <- coef(lm1)["shall"]
# Check out the SEs..
sqrt(vcovCL(lm1, cluster = ~ state_fe, type  = "HC1")["shall", "shall"])

# Initial plots... 
ggplot() + 
  geom_density(data = df, aes(x = resid, color = factor(shall))) 
ggsave(paste0(out, "1_homosked.png"), height = 3, width = 5)

filter(df, stateid == 1) %>% 
  ggplot(aes(x = year, y = resid)) + 
  geom_point()
ggsave(paste0(out, "1_state_1_time_resid.png"), height = 3, width = 5)

##############################################3
# 1.1 Try bootstrapping in all the different ways 

# Non-parametric bootstrap
npm_boot <- function(df, i, formula){
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  tibble(i = i, value = coef(reg)["shall"])
}
draws_npm <- map_dfr(seq(1:B), npm_boot, df = df, formula = reg1)
sd_npm <- sd(draws_npm$value)

# Try the residual bootstrap... 
resid_boot <- function(df, i , formula, lm){
  df$log_vio <- predict(lm) + sample(lm$residuals)
  reg        <- lm(data = df, as.formula(formula))
  tibble(i = i, value = coef(reg)["shall"])
}
draws_resid <- map_dfr(seq(1:B), resid_boot, df = df, formula = reg1, lm = lm1)
sd_resid <- sd(draws_resid$value)

# Cluster bootstrap...
cluster_boot <- function(df, i, beta_hat, formula){
  
  draw <- df %>% 
    group_nest(state_fe) %>% 
      slice_sample(prop = 1, replace = TRUE) %>% 
    unnest(cols = c(data))
  
  fit  <- lm(formula, data = draw)
  beta <- coef(fit)["shall"]
  se   <- sqrt(vcovCL(fit, cluster = draw$state_fe, type  = "HC1")["shall", "shall"])
  T    <- sqrt(nrow(draw)) * (beta - beta_hat) / se
  tibble(i = i, beta = beta, se = se, T = T)
}
draws_cluster <-  map_dfr(1:1000, cluster_boot, df = df, 
                          beta_hat = beta_hat, formula = reg1)
sd_cluster <- sd(draws_cluster$beta)

# Wild bootstrap
wild_boot <- function(df, i, formula, lm){
  
  df$omega   <- sample(c(-1,1), nrow(df) ,replace=T)
  df$log_vio <- predict(lm) + df$omega * lm$residuals
  reg        <- lm(data = df, as.formula(formula))
  tibble(i = i, value = coef(reg)["shall"])
  
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
  tibble(i = i, value = coef(reg)["shall"])
}

draws_wild_cluster <- map_dfr(seq(1:B), wild_cluster_boot, 
                              df = df, formula = reg1, lm = lm1)
sd_wild_cluster <- sd(draws_wild_cluster$value)

###########################################################
# 1.2 Plot the bootstrap outputs 

# Bind all the data together
plot_df <- draws_npm %>% 
  mutate(Bootstrap = paste0("Non-Parametric, sd = ", round(sd_npm, 4))) %>% 
  bind_rows(
    draws_resid %>% 
      mutate(Bootstrap = paste0("Residual, sd = ", round(sd_resid, 4))) 
  ) %>% 
  bind_rows(
    draws_cluster %>% select(i, value = beta) %>% 
      mutate(Bootstrap = paste0("Clustered, sd = ", round(sd_cluster, 4)))
  ) %>% 
  bind_rows(
    draws_wild %>% 
      mutate(Bootstrap = paste0("Wild, sd = ", round(sd_wild, 4)))
  ) %>% 
  bind_rows(
    draws_wild_cluster %>% 
      mutate(Bootstrap = paste0("Wild Cluster, sd = ", round(sd_wild_cluster, 4)))
  )

plot_df %>% ggplot() + 
    geom_density(aes(x = value, color = Bootstrap)) + 
    ggtitle(paste0(B, " Bootstrap Draws"))
ggsave(paste0(out, "1_bs_comparisons.png"), height= 5, width = 9)

###########################################################
# 1.3 Construct CIs using the cluster bootstrap outputs

# 1.3.1 Effron CI
# Effron... 
p_vals <- c(0.025, 0.975)
effron_ci <- quantile(draws_cluster$beta, p_vals)
effron_ci

# 1.3.1 Percentile-T CI
percentile_ci <- 
  beta_hat - quantile(draws_cluster$T, p_vals)[c(2,1)] * sd_cluster / sqrt(N)
percentile_ci

# Plot the density function, and the CIs 
ggplot() + 
  geom_density(data = draws_cluster, aes(x = beta)) + 
  geom_vline(xintercept = effron_ci, color = "red") + 
  geom_vline(xintercept = percentile_ci, color= "blue") + 
  geom_vline(xintercept = beta_hat, color  = "green") + 
  
ggsave(paste0(out, "1_CIs.png"), height = 5, width = 7)


