# Code for pset 2

############################################
# 0 Set up 
############################################
rm(list = ls())
pacman::p_load(purrr, dplyr, ggplot2, latex2exp)
theme_set(theme_minimal())
out <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/sims/exercises/2_AR_ts/")
set.seed(1)

############################################
# Problem 1
############################################

# Generate data
N <- 100
df <- 
  tibble(u = runif(N, 0, 1)
                 ) %>% 
  mutate(t = row_number(), 
         k = case_when(
                        u <= 0.25 ~  1, 
             0.25 < u & u <= 0.5  ~  2,
             0.5  < u & u <= 0.75 ~  3,
                        u  > 0.75 ~  4), 
         y = sin(0.5 * pi * (k + t))
         ) %>% 
  dplyr::select(-u) 
# Plot
ggplot(df) + 
  geom_line(aes(x = t, y = y)) + 
  geom_hline(yintercept = mean(df$y, na.rm = TRUE), 
             color = "red", alpha = 0.3)
ggsave(paste0(out, "p1_ts.png"), height = 4, width = 5)

png(paste0(out, "p1_acf.png"))
  acf(df$y, main = "ACF")
dev.off()

############################################
# Problem 2
############################################

N <- 100
df <- 
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>% 
  mutate(t = row_number(), 
         y = epsilon + 2 * lag(epsilon) + 0.64 * lag(epsilon, 2))

ggplot(df) + 
  geom_line(aes(x = t, y = y)) + 
  geom_hline(yintercept = mean(df$y, na.rm = TRUE), 
             color = "red", alpha = 0.3)
ggsave(paste0(out, "p2_ts.png"), height = 4, width = 5)

png(paste0(out, "p2_acf.png"))
  acf(df$y[!is.na(df$y)], main = "ACF")
dev.off()

############################################
# Problem 3
############################################

N <- 100
df <- 
  tibble(epsilon = rnorm(N, mean = 0, sd = 1)) %>% 
  mutate(t = row_number(), 
         y = epsilon + 1.5 * lag(epsilon) + 0.5 * lag(epsilon, 2))

ggplot(df) + 
  geom_line(aes(x = t, y = y)) + 
  geom_hline(yintercept = mean(df$y, na.rm = TRUE), 
             color = "red", alpha = 0.3)
ggsave(paste0(out, "p3_ts.png"), height = 4, width = 5)

png(paste0(out, "p3_acf.png"))
  acf(df$y[!is.na(df$y)], main = "ACF")
dev.off()




