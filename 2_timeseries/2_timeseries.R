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
# Problem 1
############################################

# Generate data
N <- 1000000
df1 <-
  tibble(u = runif(N, 0, 1)) %>%
  mutate(
    t = row_number(),
    k = case_when(
      u <= 0.25 ~ 1,
      0.25 < u & u <= 0.5 ~ 2,
      0.5 < u & u <= 0.75 ~ 3,
      u > 0.75 ~ 4
    ),
    y = sin(0.5 * pi * (k + t))
  ) %>%
  dplyr::select(-u)
# Plot
create_plots(df1, 1)

# calculate empirical autocorrelation
exp_t_t_tau <- function(df, tau){
  print(paste0("R(",tau, ") = ", mean(df$y * lag(df$y, tau), na.rm = TRUE)))
  mean(df$y * lag(df$y, tau), na.rm = TRUE)
}
exp_t_t_tau(df = df1, tau = 4)
map_dbl(1:16, exp_t_t_tau, df = df1) %>% plot()

s <- seq(1:8)
plot(s, cos(0.5*pi*s))

############################################
# Problem 2
############################################

N <- 100
df2 <-
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>%
  mutate(
    t = row_number(),
    y = epsilon + 2 * dplyr::lag(epsilon) + 0.64 * dplyr::lag(epsilon, 2)
  )
# Plot
create_plots(df2, 2)

# Compare the plots for the fundamental and non-fundamental version
df2 <- df2 %>% 
  mutate(nebla = epsilon * sqrt(2.56), 
         yn = nebla + 1.025* dplyr::lag(nebla,1) + 
           0.25 * dplyr::lag(nebla, 2))

ggplot(df2[1:20,]) + 
  geom_line(aes(x = t, y = y), color = "red") + 
  geom_line(aes(x = t, y = yn), color = "blue")

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
      var <- var + gamma(n, b[1], b[2]) * gamma(m, b[1], b[2]) * ACF[abs(n - m ) + 1]
    }
  }
  data.frame(N = N, var = var)
}

b <- c(-0.4, -0.625)
ACF <- c(5.4096, 3.28, 0.64)

plot_df2 <- map_dfr(seq(1,10), prediction_variance, ACF = ACF, b = b)
ggplot(plot_df2) + 
  geom_point(aes(x = N, y = var)) + 
  geom_hline(yintercept = 3.2, color = "red") + ylab("Forecast error variance")

ggsave(paste0(out, "p", 2, "_prediction_variance.png"), height = 4, width = 5)
for (i in 1:3) print(gamma(i, b[1], b[2]))

############################################
# Problem 3
############################################

N <- 100
df3 <-
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>%
  mutate(
    t = row_number(),
    y = epsilon + 1.5 * lag(epsilon) + 0.5 * lag(epsilon, 2)
  )
# Plot
create_plots(df3, 3)
