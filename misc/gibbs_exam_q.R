# Implements the EC0517 Exam question on gibbs sampling for censored data

library(tidyverse)

# Create a dataset
mu     <- 0.5
sample <- data.frame(Y = rnorm(n = 1000, mean = mu, sd = 1))

# Censor the data below zero
df <- 
  tibble(sample) %>% 
  mutate(y_obs = ifelse(Y<0, 0, Y))

# Plot the data, make sure its what we expect
ggplot(df) +
  geom_histogram(aes(x = y_obs))
ggplot(df) +
  geom_histogram(aes(x = Y))

gibbs <- function(df, iterations, initial_mu){
  
  # Step 0 create needed objects 
  N       <- length(df$y_obs)
  n_below <- length(df$y_obs[df$y_obs==0])
  mu      <- initial_mu
  mus     <- c(mu)
  Ys      <- data.frame(i=0, Ys=df$y_obs[df$y_obs==0])
  
  for (i in 1:iterations){
    
    # Step 1: draw Ys conditional on mu
    drawsY <- c()

    for (ii in 1:n_below){
      j <- 0
      while(j >= 0){
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
# Run the sampler
results <- gibbs(df, 10, 0.1)
mean(results$mus[2:length(results$mus)])

# Plot results, see if it makes sense... 
plot.ts(results$mus)
hist(results$mus)

results$Ys %>% 
  filter(i == 1000) %>% 
  ggplot() + 
  geom_histogram(aes(x = Ys), 
                 color = "blue", alpha = 0.5) + 
  geom_histogram(data = df, aes(x = Y), 
                 color = "red", alpha = 0.5)




