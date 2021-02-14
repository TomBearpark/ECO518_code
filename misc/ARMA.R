# slide 9 - two MA processes
pacman::p_load(tidyr, ggplot2, dplyr)
df <- 
  tibble(epsilon = rnorm(1000, mean = 0, sd = 2)) %>% 
  mutate(nebla = 0.5 * epsilon, t = row_number())

var(df$epsilon) == var(df$nebla) * 4

df <- df %>% 
  mutate(Ya = epsilon + 0.5 * lag(epsilon), 
         Yb = nebla + 2 * lag(nebla))

df %>% 
  pivot_longer(cols = c(Ya, Yb)) %>% 
  ggplot() + 
  geom_line(aes(x = t, y = value, color = name))

acf(df$Ya[!is.na(df$Ya)])
acf(df$Yb[!is.na(df$Yb)])


t <- seq(1,50)
y <- sin(0.5*pi*t)
plot.ts(t,y)

# Slide 18
