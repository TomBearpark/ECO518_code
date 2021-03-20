# testing out a parametric boostrap for normal data

library(purrr)
library(ggplot2)
library(patchwork)

boostrap <- function(N,M){
  pop <- rnorm(300000)
  
  samp <- sample(pop, N)
  
  mu <- mean(samp)
  sigma <- sd(samp)
  
  draw <- function(i, sample) tibble(
    i = i, draws = rnorm(length(sample), mean = mu, sd = sigma))
  
  df <- purrr::map_dfr(seq(1,M), draw, sample = samp)
  
  p <- df %>% 
    group_by(i) %>% 
    summarise(mean = mean(draws)) %>% 
    ggplot() + 
    geom_density(aes(x = mean))
  
  # we expect it to be distributed normal(mu, sigma_2 / N)
  q <- ggplot(data.frame(x = rnorm(10000, 0, sd = 1 / sqrt(N)))) + 
    geom_density((aes(x = x)))
  p + q
}
boostrap(1000, 20000)

# Randomization inference example from the slides
s1 <- c(rep(1, 4), 5)
s2 <- rep(0.9, 5)

sample_mean_diff <- mean(s1) - mean(s2)

s <- c(s1, s2)
get_draw_mean <- function(i, s, sample_mean_diff){
  treated_index <- sample(1:10, 10)
  s1 <- s[treated_index[1:5]]
  s2 <- s[treated_index[6:10]]
  
  data.frame(i = i, mean_diff = mean(s1) - mean(s2), 
             bigger_diff = ifelse(
    mean(s1) - mean(s2) >= sample_mean_diff, 1, 0
  ))
}
map_dfr(seq(1:5000), get_draw_mean, s = s, 
        sample_mean_diff = sample_mean_diff) %>% 
  summarise(alpha = mean(bigger_diff))

s1 <- c(rep(0.9, 4), 5)
s2 <- rep(1.1, 5)
sample_mean_diff <- mean(s1) - mean(s2)
s <- c(s1, s2)

map_dfr(seq(1:5000), get_draw_mean, s = s, 
        sample_mean_diff = sample_mean_diff) %>% 
  summarise(alpha = mean(bigger_diff))

