# ECO518 PSet 5: VARs


####################################################################
# 0. Set up
####################################################################

dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/term2/", 
        "ECO518_Metrics2/sims/exercises/5_VAR")
out <- paste0(dir, "/out/")
setwd(dir)
if(!require(VARpack2019)) 
  install.packages("VARpack2019", type = "source", repos = NULL)
library(VARpack2019)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggfortify)
library(stargazer)
library(xtable)
theme_set(theme_bw())

load("pdat.RData")
df <- pdat

####################################################################
# 1. Plot the data 
####################################################################

# AIC suggests 5 lags, btw
pacman::p_load(vars)
var_aic <- VAR(df, type = "const", lag.max = 10, ic = "AIC", season = 4)
var_aic$p

p1 <- autoplot(df, scales = "fixed")
ggsave(p1, file = paste0(out, "1_ts_faceted.png"), height = 6, width = 9)
p2 <- autoplot(df, facets = FALSE) + labs(color = "Time Series")
p2
ggsave(p2, file = paste0(out, "1_ts.png"), height = 5, width = 7)

####################################################################
# 2. Estimate a VAR
####################################################################
# Estimate a 3 lag VAR
?mgnldnsty
# Note, had to give a vprior argument to get this to work. 

VAR_and_forecast <- function(lags, df, w_prior = 0){
  
  if(lags == 3) df <- df[-(1:6), ]
  if(! lags %in% c(3,9)) stop("Not implemented for those lags")
  
  # Estimate the VAR
  var <- mgnldnsty(ydata = df, 
                    lags = lags,
                    xdata = NULL,
                    const = TRUE, 
                    vprior = list(sig = 1, w = w_prior))

  # Format initial condition
  y0 <- tail(df, lags)
  
  # Return a forecast object also
  f <- fcast(y0 = y0,         # values of time series to forecast from
             By = var$var$By,     # fitted AR coefficients from part 1
             Bx = var$var$Bx,     
             xdata = NULL,    # exogenous variables (we don't have any)
             horiz = 100,     # forecast horizon
             const = TRUE, 
             shocks = NULL)
  
  # Return a formatted data-frame for plotting 
  fdf <- f %>% 
    data.frame() %>%  
    mutate(n = row_number(), lags = paste(lags)) %>% 
    pivot_longer(cols = c("raw", "intermediate", "final")) 
  
  return(list(f = f, var=var, fdf = fdf))
}

VAR3 <- VAR_and_forecast(3, df)
VAR9 <- VAR_and_forecast(9, df)

# Compare to version with different prior
VAR3_ERIC <- VAR_and_forecast(3, df, w_prior = 20)

VAR3_ERIC$var$var$By[,,1] == VAR3$var$var$By[,,1]

VAR3$var$prior 
VAR3$var$var$Bx 

# Get the coefficients 
convert_to_latex_matrix <- function(mat, dim = 3, digits = 4){
  mat %>% 
  xtable(align = rep("", dim + 1), digits = digits) %>% 
  print(tabular.environment="bmatrix", digits = digits,
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
}
for(i in 1:3)
  convert_to_latex_matrix(VAR3$var$var$By[,,i] )

for(i in 1:9)
  convert_to_latex_matrix(VAR9$var$var$By[,,i] )

# 2.1. Plot the forecasts
plot_df <- df %>% 
  as_tibble() %>% 
  mutate(n = row_number()-dim(df)[1], lags = "historic") %>% 
  pivot_longer(cols = c("raw", "intermediate",    "final")) %>% 
  bind_rows(VAR9$fdf) %>% 
  bind_rows(VAR3$fdf)

plot_df %>% filter(n>-100) %>% 
  ggplot() + 
  geom_line(aes(x = n, y = value, color =lags)) + 
  facet_wrap(~name, scales = "fixed")
ggsave(paste0(out, "2_var_forecast_plots.png"), height = 5, width = 8)  

# 2.2. Compare the log odds ratio, print a table for latex
data.frame(VAR = c(3, 9, "Odds Ratio"), 
           Value = c(VAR3$var$w , VAR9$var$w, VAR3$var$w / VAR9$var$w)) %>% 
  stargazer(summary = FALSE)

####################################################################
# 3. Check the roots
####################################################################
T <- dim(df)[1]
?sysmat

# VAR3
mat3 <- sysmat(By=VAR3$var$var$By)
convert_to_latex_matrix(mat3, dim = 9, digits = 7)

mat3_eigen <- eigen(mat3)
mat3_eigen$values[abs(mat3_eigen$values - 1) < 1/T]

# VAR9
mat9 <- sysmat(By=VAR9$var$var$By)
mat9_eigen <- eigen(mat9)
mat9_eigen$values[abs(mat9_eigen$values - 1) < 1/T]

####################################################################
# 4. Plot the IR Plots
####################################################################
resp3 <- impulsdtrf(VAR3$var$var, nstep = 48)
plotir(resp3, file = paste0(out, "4_VAR3_IRF.pdf"), main = "IR Plot VAR 3")

resp9 <- impulsdtrf(VAR9$var$var, nstep = 48)
plotir(resp9, file = paste0(out, "4_VAR9_IRF.pdf"), main = "IR Plot VAR 9")

####################################################################
# 5. VECM
####################################################################
pacman::p_load(tsDyn)
est_VECM <- VECM(df, lag = 3, 
                  r = 1, include = "none", estim = "2OLS")
est_VECM %>% summary()
toLatex((est_VECM), parenthese="Pvalue")
coefB(est_VECM)

dfV <- 
  predict(est_VECM, n.ahead = 100) %>% 
  data.frame() %>% 
  mutate(model = "VECM", n= row_number(), lags= "VECM") %>% 
  pivot_longer(cols = c("raw", "intermediate", "final")) %>% 
  bind_rows(plot_df %>% filter(lags %in% c("historic", "3"), n > 0) ) %>% 
  mutate(model = ifelse(is.na(model), "VAR3", model))

ggplot(dfV)+
  geom_line(aes(x = n, y = value, color = name)) +
  facet_wrap(~model)
ggsave(paste0(out, "5_vecm.png"))

####################################################################
# 6. State space
####################################################################



