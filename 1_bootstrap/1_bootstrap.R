# ECO519: Exercise 1, Bootstrapping and Randomization Inference

# 0. Set up environment and paths
rm(list = ls())

root <- paste0("/Users/tombearpark/Documents/princeton/", 
               "1st_year/term2/ECO518_Metrics2/")
out <- paste0(root, "sims/exercises/1_bootstrap")


# 1. 
data <- c(3.66, 1.00, -0.87, 2.90, -0.80, 3.20, 1.69, -3.53, 3.22, 3.53)

mean <- 0
var <- 1
n <- 100

data <- rnorm(n = n, sd = var, mean = mean)

hist(data)
mean(data)

# length is 1.96 * 2

# through bootstrapping, can i get CIs right?
sample_and_mean <- function(i, data){
  data.frame(i = i, mean = sample(data, length(data), replace = TRUE))
}

# 1. Gak
draws <- map_dfr(1:100000, sample_and_mean, data = data)

draws <- draws %>% arrange(mean)

mean(draws$mean)
sd(draws$mean)
hist(draws$mean)
