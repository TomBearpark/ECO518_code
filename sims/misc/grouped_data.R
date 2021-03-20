library(haven)
library(tidyverse)

df <- read_dta("http://dss.princeton.edu/training/Panel101.dta")

# Basic OLS
ols <-lm(y ~ x1, data=df)
summary(ols)

# Fixed effects as dummies
dummy.fe <- lm(y ~ x1 + factor(country), data=df)

# Fixed effects as difference from group means
df <- 
  df %>% 
  group_by(country) %>% 
  mutate(mean.x1 = mean(x1), 
         mean.y = mean(y), 
         demeaned.x1 = x1 - mean.x1, 
         demeaned.y = y - mean.y) %>% 
  ungroup()

demeaned.fe <- lm(demeaned.y ~ demeaned.x1, data=df)

# Estimate in first differences
df <- 
  df %>% 
  group_by(country) %>% 
  mutate(fd.y = y - lag(y), 
         fd.x1 = x1 - lag(x1)) %>% 
  ungroup()

fd.lm <- lm(fd.y ~ fd.x1, data = df)
