library(tidyverse)

mu <- 1
sample <- data.frame(Y = rnorm(n = 1000, mean = mu, sd = 1))
df <- tibble(sample) %>% 
  mutate(y_obs = ifelse(Y<0, 0, Y))

ggplot(df)+
  geom_histogram(aes(x = y_obs))


gibbs <- function(df, iterations, initial_mu){
  
  # 0 get objects together
  N <- length(df$y_obs)
  n_below <- length(df$y_obs[df$y_obs==0])
  mu <- initial_mu
  mus <- c()
  Ys <- data.frame(i=0, Ys=0)
  
  for (i in 1:iterations){
    
    # Step 1: draw Ys conditional on mu
    drawsY <- c()
    for (ii in 1:n_below){
      j <- 0
      while(j>=0){
        j <- rnorm(1, mu, 1)
      }
      drawsY <- c(drawsY, j)
    }
    
    Y <- c(df$y_obs[df$y_obs != 0], drawsY)
    
    # Step 2: draw mu conditional on Y
    mu <- rnorm(1, mean = mean(Y), 1 / N)
    
    # Step 3: save results for trace plots
    mus <- c(mus, mu)
    
    Ys <- bind_rows(Ys, data.frame(i = i, Ys = Y))
  }
  return(list(mus = mus, Ys = tibble(Ys)))
}
results <- gibbs(df, 1000, 0.1)

# Plot results, see if it makes sense... 
plot.ts(results$mus)
mean(results$mus)

results$Ys %>% 
  filter(i == 100) %>% 
  ggplot() + 
  geom_density(aes(x = Ys))
