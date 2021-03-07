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
library(ggplot2)
library(ggfortify)
library(patchwork)
library(tidyr)

theme_set(theme_bw())

load("pdat.RData")
df <- pdat

####################################################################
# 1. Plot the data 
####################################################################

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

VAR_and_forecast <- function(lags, df){
  
  if(lags == 3) df <- df[-(1:6), ]
  if(! lags %in% c(3,9)) stop("Not implemented for those lags")
  
  # Estimate the VAR
  var <- mgnldnsty(ydata = df, 
                    lags = 3,
                    xdata = NULL,
                    const = TRUE, 
                    vprior = list(sig = 1, w = 0)
                    )

  # Format initial condition
  y0 <- df[(1:lags), ]
  
  # Return a forecast object also
  f <- fcast(y0 = y0,         # values of time series to forecast from
             By = var$var$By,     # fitted AR coefficients from part 1
             Bx = var$var$Bx,     
             xdata = NULL,    # exogenous variables (we don't have any)
             horiz = 100,     # forecast horizon
             const = TRUE, 
             shocks = NULL)
  return(list(f = f, var))
}
VAR_and_forecast(10, df)

VAR3 <- VAR_and_forecast(3, df)

f %>% 
  data.frame() %>%  
  mutate(n = row_number()) %>% 
  pivot_longer(cols = c("raw", "intermediate",    "final")) %>% 
  ggplot() + 
  geom_line(aes(x = n, y = value, color = name))

plot(f)
fcast(var3$var)

