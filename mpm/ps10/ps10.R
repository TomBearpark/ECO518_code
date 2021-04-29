rm(list = ls())
if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, # data wrangling
               xtable,     # printing latex tables 
               patchwork, sandwich, stargazer, lmtest)

root <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
               "ECO518_Metrics2")
dir  <- paste0(root, "/mpm/excercises/")
code <- paste0(root, "/ECO518_code/mpm")
out  <- paste0(dir, "ps10/out/");dir.create(out, showWarnings = FALSE)

setwd(dir);set.seed(1);theme_set(theme_bw())

source(paste0(code, "/ps10/ps10_funcs.R"))

###############################################################
# Question 1
###############################################################

power <- function(a, sigma1, sigma2, N1, N2){
  gamma <- a / (sqrt(sigma1^2) / N1 + sqrt(sigma2^2) / N2)
  z <- qnorm(.975)
  pnorm(gamma - z) + pnorm(-gamma - z)
}
a <- seq(-3,3,0.1)
tibble(a = a, power = power(a, 2,2,10,10)) %>% 
  ggplot() + geom_line(aes(x = a, y = power)) + 
  geom_hline(yintercept = 0, color = "red", alpha = 0.4) + 
  ggsave(paste0(out, "1_power_function_example.png"), height = 4, width = 6)

power_df_N <- function(N, sigma = 2){
  N1 <- N2 <- N / 2
  sigma1 <- sigma2 <- sigma
  a <- seq(-5,5,0.1)
  tibble(N = N, a = a, power = power(a, sigma1,sigma2, N1, N2))
}
power_df_sigma <- function(sigma){
  N1 <- N2 <- 20 / 2
  sigma1 <- sigma2 <- sigma
  a <- seq(-5,5,0.1)
  tibble(sigma = sigma, a = a, power = power(a, sigma1,sigma2, N1, N2))
}
# cray graph for fun... 
(map_dfr(1:100, power_df_N) %>% 
ggplot(aes(group = N)) + 
  geom_line(aes(x = a, y = power, color = N)) + 
  scale_color_viridis_c()) + 
(map_dfr(1:100, power_df_sigma) %>% 
   ggplot(aes(group = sigma)) + 
   geom_line(aes(x = a, y = power, color = sigma)) + 
   scale_color_viridis_c()) + 
  ggsave(paste0(out, "1_power_function_N.png"), height = 4, width = 9)

###############################################################
# Question 2
###############################################################
df <- read_csv(paste0(dir, "/ps10/jtpa.csv"))

# i
lm1 <- lm(data = df, earnings ~ treatment)
sd1 <- sqrt(diag(vcovHC(lm1, type = "HC1")))
coeftest(lm1, vcov = vcovHC(lm1, type = "HC1")) %>% stargazer()

# ii
df$age <- case_when(
  df$age2225 == 1 ~ "22 25",
  df$age2629 == 1 ~ "26 29", 
  df$age3035 == 1 ~ "30 35", 
  df$age3644 == 1 ~ "36 44", 
  df$age4554 == 1 ~ "45 54")
df$age <- ifelse(is.na(df$age), "55 78", df$age)

lm2 <- lm(data = df, earnings ~ treatment + factor(age))
coeftest(lm2, vcov = vcovHC(lm2, type = "HC1")) %>% stargazer()
sd2 <- sqrt(diag(vcovHC(lm2, type = "HC1")))

# Print for latex
stargazer(lm1, lm2, se = list(sd1, sd2))

# iii
power <- function(N, df){
  z <- qnorm(.975)
  a <- 1000
  sigma1 <- df$earnings[df$treatment == 1] %>% sd()
  sigma0 <- df$earnings[df$treatment == 0] %>% sd()
  N1 <- 2/3 * N
  N2 <- 1/3 * N
  gamma <- a / (sqrt(sigma1^2 / N1 + sigma0^2 / N2))
  pnorm(gamma - z) + pnorm(-gamma - z)
}
res <- map_dbl(1:10000, power, df)

plot_df <- tibble(N = 1:10000, power = res)
plot_df %>% ggplot() + 
  geom_line(aes( x = N, y = power )) + 
  geom_hline(yintercept = 0.8, color = "red") + 
  geom_hline(yintercept = 0, color = "black") + 
  ggsave(file = paste0(out, "2_power_N.png"), height = 5, width = 6)

plot_df$N[plot_df$power == min(res[res>.80])]

###############################################################
# Question 3
###############################################################

df <- 
  read_csv(paste0(dir, "ps10/draft.csv")) %>% 
  mutate(c = 1)

df %>% ggplot() + geom_density(aes(x = lwage, color = factor(yob)))
df %>% ggplot() + geom_density(aes(x = lwage, color = factor(veteran))) + 
  facet_wrap(~yob)

get_2sls <- function(yob, df){
  
  df <- filter(df, yob == !!yob)
  Y <- as.matrix(df[c("lwage")])
  X <- as.matrix(df[c("c", "veteran")])
  R <- as.matrix(df[c("c", "draftelig")])
  
  res <- estimate_2SLS(Y, X, R) %>% mutate(yob = yob)
  res$parameter <- c("constant", "veteran")
  res
}
yobs <- unique(df$yob)
results <- map_dfr(yobs, get_2sls, df = df) %>% 
  pivot_longer(cols = c(coefficient, sd)) %>% 
  pivot_wider(id_cols = c(parameter, name), 
                        names_from = yob, values_from = c(value), 
              names_prefix = "yob_") 

results %>% xtable() %>% print()

# Fraction of compliers, never takers, always takers 
fracs <- function(df, yob){
  df <- filter(df, yob == !!yob)
  N <- dim(df)[1]
  N1 <- 1*(df$draftelig == 1) %>% sum()
  N0 <- 1*(df$draftelig == 0) %>% sum()
  p11 <- sum(df$veteran[df$draftelig == 1]) / N1
  p10 <- sum(df$veteran[df$draftelig == 0]) / N0
  tibble(
    yob = yob, N = N, 
    compliers = p11-p10, 
    always_take = p10, 
    never_take = 1 - p11
  )
}
fracs_df <- map_dfr(yobs, fracs, df = df )
fracs_df %>% xtable() %>% print(include.rownames=FALSE)
