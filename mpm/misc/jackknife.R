jackknife <- function(func, data)
{
  if(length(formalArgs(mean)[formalArgs(mean) != "..."]) != 1 | !is.vector(data)
    ) stop("This function is only implemented for a univariate situation")
  
  N <- length(data)
  
  theta_minus_is <- c()
  for(i in 1:N){
    df <- data[-i]
    theta_minus_i  <- func(df)  
    theta_minus_is <- c(theta_minus_i, theta_minus_is)
  }
  
  theta_hat_bar <- mean(theta_minus_is)
  deviation     <- theta_minus_is - theta_hat_bar
  sqrt(((N-1) / N) * sum(deviation^2))
}
