library(dplyr)
library(ggplot2)

mat <- matrix(
  c(1.3, -0.5, 
    1,0), nrow = 2
)
invert <- function(x) 1 / x
eigen(mat)$values %>% invert

# Simulate the series 
df <- data.frame(epsilon = rnorm(1000, 0, 2), 
                 X = 0)
df$X[1] <- 10 + df$epsilon[1] 
df$X[2] <- 2 + 1.3 * df$X[1] +  df$epsilon[2] 

for(i in 3:length(df$X)){
  df$X[i] <- 2 + 1.3 * df$X[i-1]- 0.5 * df$X[i-1] +  df$epsilon[i] 
}
plot.ts(df$X)
mean(df$X)

# Calculate the variance: 
# 1. using kronecker product
B <- mat
I <- diag(1, dim(kronecker(B,B)))
sigma <- matrix(c(4,0,0,0), nrow = 2)
sigmaBar <- matrix(sigma, nrow = 4)
Ry0Bar <- solve(I - kronecker(B,B)) %*% sigmaBar


# 2. using doubling algorithm
omega <- sigma
W <- B

for (i in 1:10){
  omega <- omega + W %*% omega %*% t(W)
  W <- W %*% W
  print(omega[1,1])
}
omega


mat <- matrix(
  c(1.1,0,-0.4,-0.2, 
    0.3, 1.4, 0, -0.2, 
    1,0,0,0,
    0,1,0,0), nrow= 4, byrow = TRUE
)
mat
decomposition <- eigen(mat)
decomposition$vectors %>% solve



epsilon <- rnorm(100)
epsilon <- c(0,1,rep(0, 20))
phi1 <- 1
phi2 <- 2
y <- c()
for (i in 2:10){
  y[i] <- phi1*epsilon[i] + phi2*epsilon[i-1]
}
plot.ts(y)

mat <- matrix(c(3,2,2, 3), nrow = 2)
mat <- matrix(c(2,3,3, 2), nrow = 2)
mat
t(chol(mat)) %*% chol(mat)
eigen(mat)$vectors 
