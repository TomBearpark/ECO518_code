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
