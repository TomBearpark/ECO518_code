# Code for PSET 2
# Creates visualizations of the time series and ACFs for each of the
# models.

############################################
# 0 Set up
############################################
rm(list = ls())
pacman::p_load(purrr, dplyr, ggplot2, latex2exp)
theme_set(theme_minimal())
out <- paste0(
  "/Users/tombearpark/Documents/princeton/1st_year/term2/",
  "ECO518_Metrics2/sims/exercises/2_AR_ts/"
)
set.seed(1)

# Function for saving a time series and acf for each model
create_plots <- function(df, problem) {
  ggplot(df) +
    geom_line(aes(x = t, y = y)) +
    geom_hline(
      yintercept = mean(df$y, na.rm = TRUE),
      color = "red", alpha = 0.3
    )
  ggsave(paste0(out, "p", problem, "_ts.png"), height = 4, width = 5)

  png(paste0(out, "p", problem, "_acf.png"))
  acf(df$y, main = "ACF")
  dev.off()
}

############################################
# Problem 2
############################################

# Visualise the model 
N <- 100
df2 <-
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>%
  mutate(
    t = row_number(),
    y = epsilon + 2 * dplyr::lag(epsilon) + 0.64 * dplyr::lag(epsilon, 2)
  )
# Plot
create_plots(df2, 2)

# Function to calculate the coefficient in the lag 
# polynomial 
gamma <- function(n, b1, b2){
  coef <- 0
  for (i in 0:n){
    coef <- coef + (b1)^(n - i)*(b2)^(i)
  }
  coef
}

# Function to calculate prediction variance, see overleaf for details 
prediction_variance <- function(N, ACF, b){
  ACF <- append(ACF, rep(0, N))
  var <- ACF[1]
  for (n in 1:N) {
    var <- var + 2 * gamma(n, b[1], b[2]) * ACF[n + 1]
    for (m in 1:N){
      var <- var + 
        gamma(n, b[1], b[2]) * gamma(m, b[1], b[2]) * ACF[abs(n - m ) + 1]
    }
  }
  data.frame(N = N, var = var)
}

# Calculate the variance of each estimator, plot as a function of N
b <- c(-0.4, -0.625)
ACF <- c(5.4096, 3.28, 0.64)
plot_df2 <- map_dfr(seq(1,10), prediction_variance, ACF = ACF, b = b)
ggplot(plot_df2) + 
  geom_point(aes(x = N, y = var)) + 
  geom_hline(yintercept = 3.2, color = "red") + ylab("Forecast error variance")
ggsave(paste0(out, "p", 2, "_prediction_variance.png"), height = 4, width = 5)

# Print out the coefficients of the selected model to copy into overleaf
for (i in 1:3) print(gamma(i, b[1], b[2]))

# check the results of the prediction - get predicted values based on our 
# model
get_prediction <- function(order, df){
  df[paste0('pred_N')] <- 0
  for (i in 1:order){
    df[paste0('pred_N')] <- df[paste0('pred_N')] - 
      gamma(i, b[1], b[2]) * dplyr::lag(df$y, i) 
  }
  data.frame(order = order, pred = df[paste0('pred_N')], y = df$y) %>% 
    mutate(error = pred_N - y)
}

# Scatter predicted vs actual values for differing lag lengths 
plot_df2 <- map_dfr(1:10, get_prediction, df = df2)
ggplot(plot_df2 %>% filter(order < 5)) +
  geom_point(aes(x = y, y = pred_N)) + 
  geom_smooth(aes(x = y, y = pred_N), alpha = 0.1)+
  facet_wrap(~order)
ggsave(paste0(out, "p", 2, "_empirical_prediction_scatter.png"), 
       height = 4, width = 5)

# Make a table of variances to copy into overleaf
plot_df2 %>% 
  group_by(order) %>% 
  summarise(var = sd(error, na.rm = TRUE)^2) %>% 
  xtable()

############################################
# Problem 3
############################################
# Generate time series from the model
N <- 100
df3 <-
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>%
  mutate(
    t = row_number(),
    y = epsilon + 1.5 * lag(epsilon) + 0.5 * lag(epsilon, 2)
  )
# Plot
create_plots(df3, 3)

# Function for getting the coefficient on the lag
gammaK <- function(n, K, b1){
  coef <- 0
  for (i in 0:n){
    coef <- coef + (-(1 - n/K))^(n - i)*(b1)^(i)
  }
  coef
}
# Find the variance for a given model
prediction_variance2 <- function(N, ACF, b){
  ACF <- append(ACF, rep(0, N))
  var <- ACF[1]
  for (n in 1:N) {
    var <- var + 2 * gammaK(n, K = N, b) * ACF[n + 1]
    for (m in 1:N){
      var <- var + 
        gammaK(n, K = N, b)  * gammaK(m, K = N, b)  * ACF[abs(n - m) + 1]
    }
  }
  data.frame(N = N, var = var)
}
# PLot prediction variance as a funciton of N
ACF <- c(3.5, 2.25, 0.5)
b <- -0.5
plot_df3 <- map_dfr(seq(1,30), prediction_variance2, ACF = ACF, b = b)
ggplot(plot_df3) + 
  geom_point(aes(x = N, y = var)) + 
  geom_hline(yintercept = 1.25, color = "red")
ggsave(paste0(out, "p", 3, "_prediction_variance.png"), 
       height = 4, width = 5)

# print coefficients to copy into latex
l <- ""
for(n in 1:12) 
  l <- paste0(l, " + ",round(gammaK(n, 12, -0.5), 4), "y_{t-", n, "}")
l

