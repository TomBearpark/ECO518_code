# Code for ECO518 MPM Pset 1 

###########################################################
library(tidyverse)
library(readxl)
theme_set(theme_bw())

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps6/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(123)

###########################################################
# Question 1
B <- 1000

df <- read_xlsx("Guns.xlsx")
df <- mutate(df, 
             log_vio = log(vio), 
             state_fe = factor(stateid))

# linear FE regression..
reg1 <- paste0("log_vio ~ shall + incarc_rate + density + avginc + pop ", 
               "+ pb1064 + pw1064 + pm1029 + state_fe")

lm1 <- lm(data = df, formula = as.formula(reg1))
df$resid <- lm1$residuals
summary(lm1)
plot(lm1)
ggplot(data = df, aes(y = resid, x = shall)) + 
  geom_point(col = 'blue') + 
  geom_abline(slope = 0)
# seems pretty homoskedastic

# Bootstrap the standard errors...

# Non-parametric bootstrap
npm_boot <- function(df, i, formula){
  n    <- nrow(df)
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}
draws_npm <- map_dfr(seq(1:B), npm_boot, df = df, formula = reg1)

# Find the sd:
sd(draws_npm$value)

# Try the residual bootstrap... 
resid_boot <- function(df, i , formula, lm){
  df$log_vio <- predict(lm) + sample(lm$residuals)
  reg <- lm(data = df, as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}
draws_resid <- map_dfr(seq(1:B), resid_boot, 
                       df = df, formula = reg1, lm = lm1)
sd(draws_resid$value)

# Cluster bootstrap...
cluster_boot <- function(df, i, formula){
  
  draw <- df %>% 
    group_by(state_fe) %>% 
      slice_sample(prop = 1, replace = TRUE) %>% 
    ungroup()
  
  reg  <- lm(data = draw, formula = as.formula(formula))
  data.frame(i = i, value = coef(reg)["shall"])
}

draws_cluster <- map_dfr(seq(1:B), cluster_boot, 
                       df = df, formula = reg1)
sd(draws_cluster$value)

# very similar - homo-skedasticity might be a good assumption 
ggplot() + 
  geom_density(data = draws_npm, aes(x = value), color = "blue") + 
  geom_density(data = draws_resid, aes(x = value), color = "red") + 
  geom_density(data = draws_cluster, aes(x = value), color = "green")




