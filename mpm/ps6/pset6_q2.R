## UNFINISHED - and not requried for the pset! 

# Question 2
N <- 1000
epsilon <- rnorm(N, 0, 1000)


# Create the underling data
X <- sample(seq(1,100), size= N, replace = TRUE)
beta <- 10000
Z_star <- 100000 + beta * X + epsilon
hist(Z_star)


# Censor the data 
Z <- ifelse(Z_star< 10^6, Z_star, NA)
hist(Z)
data <- data.frame(X = X, Z = Z)


# Step 1 - Compute the MLE estimator 

L <- function(param, data, N){
  
  gamma   <- param[1]
  sigma_2 <- param[2]
  
  L <- 0
  
  for(i in 1:N){
    
    x <- data$X[i]
    z <- data$Z[i]
    
    # print((1 / sqrt(sigma_2 * 2 * pi)))
    # * 
    # q <- exp(-1 /(2 * sigma_2) * (z - gamma * x)^2)
    # print(log(q))
    # print() 
    # 
    if(!is.na(z)){
      # browser()
      L <- L + log(1 / sqrt(sigma_2 * 2 * pi)) + 
                -1 /(2 * sigma_2) * (z - gamma * x)^2
    }else{
      L <- L + log(1 - pnorm((10^6 - gamma * x) / sqrt(sigma_2)))
    }
    # print(L)
  }
  L
}

L(param = c(10000,1000), data, N)

optim(par = c(1,1), fn = L, data = data, N = N)
