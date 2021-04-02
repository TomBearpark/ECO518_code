# code for MPM pset 7

rm(list = ls())
library(tidyverse) 
library(furrr)
theme_set(theme_bw())
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
              "ECO518_Metrics2/mpm/excercises/ps7/")
out <- paste0(dir, "out/")
setwd(dir)
set.seed(1)
# plan(multisession, workers = 8)
###########################################################
# Question 1
###########################################################

# 0 Load data, check it out
df <- read_csv(paste0(dir, "engel.csv")) %>% 
  mutate(log_inc = log(income))

ggplot(df) + 
  geom_point(aes(x = foodexp, y = log_inc))

###########################################################
# 1.(i) 
# Using a normal kernel and the normal reference rule, estimate 
# the density of log income at each value of X, where X is the 
# log income data. Plot the estimated density of log income.

# Normal kernel function
K <- function(Xi, X0, h){
  z <- (Xi - X0) / h
  (2*pi) ^ (-.5) * exp(- (1 / 2) * z^2 ) 
}

# Density function, returns the density at X0, given data X and bandwidth h
f <- function(X0, X, h){
  N <- length(X)
  fhat <- 0
  for(i in 1:N){
    fhat <- fhat + K(Xi = X[i], X0 = X0, h = h)
  }
  data.frame(X = X0, density = (1 / N*h) * fhat)
}

# Calculate h using normal reference rule 
h_star <- function(X){
  N <- length(X)
  sigma <- sd(X)
  1.059 * sigma  / N ^ (0.2)
}

X <- df$log_inc

h <- h_star(X)
density <- map_dfr(X, f, X = X, h = h)
ggplot(density) + 
  geom_point(aes( x = X, y = density)) + 
  ggtitle(paste0("bandwidth  = ", h))


###########################################################
# 1 (ii)
df <- read_csv(paste0(dir, "engel.csv")) %>% 
  mutate(log_inc = log(income))
df$share_food <- df$foodexp / df$income

Y <- df$share_food
X <- df$log_inc

m_nw <- function(X0, X, Y, h){
  K <- function(Xi) dnorm((Xi - X0) / h)
  w <- sapply(X, K) %>% rbind()
  denom <- rowSums(w)
  w <- w / denom
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = K(X0)/denom, h = h)
}

CV <- function(hvals, X, Y, mfunc = m_nw){
  CV_out <- data.frame(h = hvals, CV = 0)
  N <- length(Y)
  i <- 1
  for(h in hvals){
    m <- map_dfr(X, mfunc, X, Y, h = h)
    CV_out$CV[i] <- sum(((Y - m$mhat) / (1 - m$wi_Xi))^2) / N
    i <- 1 + i
  }
  CV_out
}

findMin <- function(CV_results) {
  CV_results$h[CV_results$CV == min(CV_results$CV)]
}

CV_nw_results <- CV(hvals = seq(0.1, 1, 0.01) , X = X, Y = Y, mfunc = m_nw )
CV_nw_results %>% plot()
CV_h_nw <- findMin(CV_nw_results)

rangeX <- seq(min(X), max(X), length.out = 100)
mhat_nw <- m_nw(rangeX, X, Y, CV_h_nw)
ggplot() + 
  geom_line(data = mhat_nw, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)

# Try and do a cool plot
plot_df <- map_dfr(c(seq(0.1, 1, 0.05), CV_h_nw), 
                   m_nw, X0 = rangeX, X = X, Y = Y)

ggplot(plot_df,  aes(x = X0, y = mhat, group = h)) + 
  geom_line(aes(color = h), alpha = 0.6)  + 
  geom_line(data = filter(plot_df, h == CV_h_nw), 
            aes(x = X0, y = mhat), color = "red") + 
  scale_color_viridis_c() + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)

###########################################################
# 1 (iii)
m_ll <- function(X0, X, Y, h){
  
  K <- function(Xi) dnorm((Xi - X0) / h)
  z <- matrix(c(rep(1, length(X0)), X0), nrow = 2)
  mid <- 0
  for(j in 1:length(X)){
    Zj  <- matrix(c(1, X[j]))
    mid <- mid + K(X[j]) * Zj %*% t(Zj)
  }
  mid <- solve(mid)
  
  weight <- function(Xi) t(z) %*% mid %*% (K(Xi) * matrix(c(1, Xi)))
  
  w <- sapply(X, weight) %>% rbind()
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = weight(X0), h = h)
}

CV_ll_results <- CV(hvals = seq(1,2, 0.01), X = X, Y = Y, mfunc = m_ll) 
CV_ll_results %>% plot()
CV_h_ll <- findMin(CV_ll_results)
mhat_ll <- map_dfr(X, m_ll, X, Y, h = CV_h_ll)

ggplot() + 
  geom_line(data = mhat_ll, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)

# Cool plot to show how the bandwidth affects the estimates 
plot_df2 <- data.frame(X0 = 0, mhat = 0, wi_Xi = 0, h = 0)
for(h in c(seq(0.5, 3, 0.1), CV_h_ll)){
  dfh <- map_dfr(rangeX, m_ll, X, Y, h = h)
  plot_df2 <- bind_rows(plot_df2, dfh)
}
plot_df2 <- tibble(plot_df2)
plot_df2 <- plot_df2[-1,] 
ggplot(plot_df2,  aes(x = X0, y = mhat, group = h)) + 
  geom_line(aes(color = h))  +
  geom_line(data = filter(plot_df2, h == CV_h_ll), aes(x = X0, y = mhat), color = "red") + 
  scale_color_viridis_c() + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2) + 
  ggtitle("Local Linear Regression Estimates")



###########################################################
# 1 (iv)

z <- function(x, p){
  z <- c(1)
  for(j in 1:p) z[j+1] <- x^j
  matrix(z)
}

m_poly <- function(X0, X, Y, h){
  
  p <- h
  mid <- 0
  for(j in 1:length(X)){
    Zj  <- z(X[j], p)
    mid <- mid + Zj %*% t(Zj)
  }
  mid <- solve(mid)
  
  weight <- function(Xi) t(z(X0, p)) %*% mid %*%  z(Xi, p)
  
  w <- sapply(X, weight) %>% rbind()
  mhat <- drop(w %*% Y)
  data.frame(X0 = X0, mhat = mhat, wi_Xi = weight(X0))
}

pvals <- seq(1,10)
CV_results_poly <- data.frame(h = pvals, CV = 0)
for(p in pvals){
  print(p)
  try(CV_results_poly$CV[p] <- CV(hvals = p, X, Y, mfunc = m_poly)$CV)
}
plot(CV_results_poly$CV)
mhat_poly <- map_dfr(rangeX, m_poly, X, Y, h = 2)
ggplot() + 
  geom_line(data = mhat_poly, aes(x = X0, y = mhat)) + 
  geom_point(data = df, aes(x = log_inc, y = share_food), alpha = .2)


###########################################################
# Question 2
###########################################################

# i)

df2 <- read_csv("nls.csv") %>% 
  rename(Y = luwe, X = educ, R = exper) %>% 
  mutate(R2 = R^2)
parametric <-  "Y ~ X + R + R2"
reg1 <- lm(df2, formula = parametric) 
summary(reg1)
coef(reg1)["X"]

# ii)

Y <- df2$Y
R <- df2$R
X <- df2$X

rangeR <- seq(min(R), max(R), length.out = 100)
hvals <- seq(0.1, 1, 0.1)

# E[X | R]
CV_X.R <- CV(hvals = seq(0.1, 1, 0.05) , X = R, Y = X, mfunc = m_nw)
CV_X.R <- filter(CV_X.R, !is.na(CV))
CV_X.R %>% plot()
CV_h_X.R <- findMin(CV_X.R)
Xhat.R <- m_nw(rangeR, X = R, Y = X, CV_h_X.R)
ggplot() +
  geom_line(data = Xhat.R, aes(x = X0, y = mhat)) +
  geom_point(data = df2, aes(x = R, y = X), alpha = .2) + 
  xlab("R") + ylab("X")

# E[Y | R]
CV_Y.R <- CV(hvals = seq(0.1, 1, 0.05) , X = R, Y = Y, mfunc = m_nw)
CV_Y.R <- filter(CV_Y.R, !is.na(CV))
CV_Y.R %>% plot()
CV_h_Y.R <- findMin(CV_Y.R)
Yhat.R <- m_nw(rangeR, X = R, Y = Y, CV_h_Y.R)
ggplot() +
  geom_line(data = Xhat.R, aes(x = X0, y = mhat)) +
  geom_point(data = df2, aes(x = R, y = X), alpha = .2) + 
  xlab("R") + ylab("Y")


# Calculate the expected values: 
Xhat.R <- m_nw(R, X = R, Y = X, CV_h_X.R)
df2$Xhat <- Xhat.R$mhat
Yhat.R <- m_nw(R, X = R, Y = Y, CV_h_Y.R)
df2$Yhat <- Yhat.R$mhat

df2 <- df2 %>% 
  mutate(y_res = Y - Yhat, 
         x_res = X - Xhat)
ggplot(df2) + 
  geom_point(aes(x = x_res, y = y_res))

reg <- lm(data = df2, formula = "y_res ~ 0 + x_res")


# Bootstrap the estimators

# Simple npm bootstrap for (i) estimator 
B <- 1000
npm_boot <- function(df, i, formula){
  draw <- slice_sample(df, prop = 1, replace = TRUE)
  reg  <- lm(data = draw, formula = as.formula(formula))
  tibble(i = i, value = coef(reg)["X"])
}
draws_npm <- map_dfr(seq(1:B), npm_boot, df = df2, formula = parametric)
sd_npm <- sd(draws_npm$value)
ggplot(draws_npm) + 
  geom_density(aes(x = value))

# NPM estimator for (ii)
boot <- function(df, i, hX, hY){
  df <- slice_sample(df, prop = 1, replace = TRUE)
  X <- df$X
  Y <- df$Y
  R <- df$R
  Xhat.R <- m_nw(R, X = R, Y = X, hX)
  df$Xhat <- Xhat.R$mhat
  Yhat.R <- m_nw(R, X = R, Y = Y, hY)
  df$Yhat <- Yhat.R$mhat
  df <- df %>% 
    mutate(y_res = Y - Yhat, 
           x_res = X - Xhat)
  reg <- lm(data = df, formula = "y_res ~ 0 + x_res") 
  tibble(i = i, value = coef(reg)["x_res"])
}
draws_boot <- map_dfr(seq(1:B), boot, df = df2, hX = CV_h_X.R, hY = CV_h_Y.R)
sd(draws_boot$value)
mean(draws_boot$value)
ggplot() + 
  geom_density(data = draws_boot, aes(x = value), color = "blue") + 
  geom_density(data = draws_npm, aes(x = value), color = "red")  
