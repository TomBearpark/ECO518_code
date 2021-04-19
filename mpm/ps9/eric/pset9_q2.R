library(abind)
library(tidyverse)
library(janitor)
library(furrr)

# Read in the data, clean up the names
df <-  read_csv("ps9/heating.csv") %>% 
  janitor::clean_names() %>% 
  mutate(i = row_number())
# Note - i use indexes 1:5 rather than 0:4 for convenience
num_choices <- 5
N <- dim(df)[1]

names(df) <- 
  c(paste0("ic_", seq(1,num_choices)), 
    paste0("oc_", seq(1,num_choices)), "choice", "i")

# Create choices matrix 
for(j in 0:4){
  df[paste0("Y_i", j+1)] <- ifelse(df$choice == j, 1, 0)
}
Y <- df[c(paste0("Y_i", 1:5))] %>% as.matrix()

# Create covariate array - N x k x NumChoices 
covariates <- c("ic", "oc")
k <- length(covariates)
get_x <- function(df, N, k, num_choices){
  X <- array(dim = c(N, k, num_choices))
  for(i in 1:5){
    X[,,i] <- as.matrix(df[c(paste0("ic_",i), paste0("oc_", i))])
  }
  X
}
X <- get_x(df, N,k, num_choices)

# Add columns for the constants - don't give one to the first one, which is 
# normalised to zero 
nC <- num_choices - 1
X <- abind(array(0, replace(dim(X), k, nC)), X, along = 2)
for (j in 1:nC) X[, j, j+1] <- 1

# Returns negative of the log likelihood
loglikCondProbit <- function(theta, X, Y, constant = TRUE){
  
  if(!constant) X <- get_x(df, N, k, num_choices)
  
  N <- dim(Y)[1]
  log_lik <- 0
  for(i in 1:N){
    X_i <- X[i,,]
    Y_i <- Y[i,]
    log_lik <- log_lik + 
      t(theta) %*% X_i %*% Y_i - log(sum(exp(t(theta) %*% X_i)))
  }
  -drop(log_lik)
}

# Run MLE to get an estimate of the paramter 
theta <- matrix(rep(0, 6))
MLE <- optim(par = theta, loglikCondProbit, X = X, Y = Y)
theta <- MLE$par
theta
# Get the standard errors using the information matrix estimator 
Omega <- matrix(0, ncol = length(theta), nrow = length(theta))
for (i in 1:N ){
  X_i <- X[i,,]
  Y_i <- Y[i,]
  p_i <- exp(t(theta) %*% X_i)/sum(exp(t(theta) %*% X_i)) 
  alpha_i <- Y_i[2:5] - p_i[2:5]
  beta_i <- (Y_i - p_i) %*% t(X_i[(nC+1):(nC + k), ])
  Delta_f <- matrix(c(alpha_i, beta_i))
  Omega <- Omega + Delta_f %*% t(Delta_f)
}
Omega <- Omega/N 
V     <- solve(Omega)
SE    <- diag(sqrt(V/N))

# iv) likelihood ratio
MLE_NC <- optim(par = c(0,0), loglikCondProbit, X = X, Y = Y, constant = FALSE)
LR <- -2*(-MLE_NC$value - (-MLE$value))

# v) Market Shares
average_p <- function(X, Y, N, theta){
  p_i <- 0
  for (i in 1:N){
    X_i <- X[i,,]
    Y_i <- Y[i,]
    p_i <- p_i + exp(t(theta) %*% X_i)/sum(exp(t(theta) %*% X_i)) 
  }
  p_i / N
}

E_Pij            <- average_p(X, Y, N, MLE$par)
X_altered        <- X
X_altered[, 6, c(2,4)] <- X[, 6, c(2,4)] * 1.1

E_Pij_altered    <- average_p(X_altered, Y, N, MLE$par)
E_Pij_altered - E_Pij

boot <- function(i, initial, X, Y, N){
  print(i)
  index <- sample(1:N, size = N, replace = TRUE)
  X <- X[index, , ]
  Y <- Y[index, ]
  MLE <- optim(par = initial, loglikCondProbit, X = X, Y = Y)
  E_Pij <- average_p(X, Y, N, MLE$par)
  X_altered <- X
  X_altered[, 6, c(2,4)] <- X[, 6, c(2,4)] * 1.1
  E_Pij_altered <- average_p(X_altered, Y, N, MLE$par)
  E_Pij_altered - E_Pij
}
boot(1, MLE$par, X, Y, N)
plan(multisession(workers = 7))
draws <- future_map(seq(1,1000), boot, initial = MLE$par, X = X, Y= Y, N = N)

draws <- do.call(rbind.data.frame, draws)
draws %>%  summarize_all(sd)
