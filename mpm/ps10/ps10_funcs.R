estimate_2SLS <- function(Y, X, R, return.omega = FALSE){
  
  N <- length(Y)
  S_rx <- t(R) %*% X / N
  S_rr <- t(R) %*% R / N
  S_ry <- t(R) %*% Y / N
  
  # calculate beta
  beta <- solve(t(S_rx) %*% solve(S_rr) %*% S_rx)  %*% 
    t(S_rx) %*% solve(S_rr) %*% S_ry
  
  # calculate SE
  u_hat <- Y - X %*% beta
  Omega <- diag(N)
  diag(Omega) <- u_hat^2
  
  G <- -1/N * t(R) %*% X
  W <-  1/N * solve(t(R) %*% R)
  O <-  1/N * t(R) %*% Omega %*% R
  
  gmm_var <-  solve(t(G) %*% W %*% G) %*% 
    t(G) %*% W %*% O %*% W %*% G %*% 
    solve(t(G) %*% W %*% G)
  
  sd_gmm <- sqrt(diag(gmm_var / N))
  if(return.omega){
    list(results = tibble(coefficient = drop(beta), sd = sd_gmm), 
         omega = O)
  }else{
    tibble(coefficient = drop(beta), sd = sd_gmm)
  }
}
