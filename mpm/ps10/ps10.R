rm(list = ls())
if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, # data wrangling
               xtable     # printing latex tables 
)
root <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
               "ECO518_Metrics2")
dir <- paste0(root, "/mpm/excercises/")
code <- paste0(root, "/ECO518_code/mpm")
out <- paste0(dir, "ps10/out/");dir.create(out, showWarnings = FALSE)

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


###############################################################
# Question 2
###############################################################
df <- read_csv(paste0(dir, "/ps10/jtpa.csv"))

df$age <- case_when(
  df$age2225 == 1 ~ "22 25",
  df$age2629 == 1 ~ "26 29", 
  df$age3035 == 1 ~ "30 35", 
  df$age3644 == 1 ~ "36 44", 
  df$age4554 == 1 ~ "45 54", 
  is.na(df$age) == TRUE ~ "55 78"
)

ggplot(data = df) + 
  geom_density(aes(x = earnings, color= factor(treatment))) + 
  facet_wrap(~age)

# regii <- lm(data = df, 
            # earnings ~ )

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
  res$parameter = c("constant", "veteran")
}
yobs <- unique(df$yob)
results <- map_dfr(yobs, get_2sls, df = df)
results %>% 
  pivot_longer(cols = c(coefficient, sd)) %>% 
  pivot_wider(id_cols = c(parameter, name), 
                        names_from = yob, values_from = c(value), 
              names_prefix = "yob_") %>% xtable() %>% print()

pacman::p_load(ivpack,ivreg)
ivreg(lwage ~ veteran | draftelig, data = filter(df, yob == 50)) %>% robust.se()

