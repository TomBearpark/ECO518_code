# ECO518 PS4

##############################################
# 0 Set up
##############################################
# Load packages
if(!require(pacman)) install.packages("pacman")
pacman::p_load(ggplot2, dplyr, sf, tigris, viridis, patchwork, sandwich, nlme)
theme_set(theme_bw())
# Set paths
dir <- paste0("/Users/tombearpark/Documents/princeton/1st_year/", 
              "term2/ECO518_Metrics2/sims/exercises/4_grouped_data/")
out <- paste0(dir, "out/")
setwd(dir)
# Load in the data
load("caschool.RData")
df <- tibble(caschool)
# load a shapefile for maps
cal <- counties(state = "California", cb = TRUE)

##############################################
# 1. Data vis
##############################################

plot_df <- 
  left_join(cal, df %>% group_by(county) %>% tally(), 
      by = c("NAME" = "county")) %>%  
  mutate(obs = ifelse(is.na(n), 0, n))

# Make sure the merge worked properly 
stopifnot(
  dim(anti_join(df %>% group_by(county) %>% tally(), 
                cal, by = c("county" = "NAME")))[1] == 0)

wrap_elements(
  (ggplot(plot_df) + geom_sf(aes(fill = n)) + 
    scale_fill_viridis(na.value = "white")) + 
  (ggplot(plot_df) + geom_histogram(aes(x = n)))) 

##############################################
# 2. Exercise 
##############################################
# Estimate a linear regression of the average test score (testscr) on student-teacher
# ratio, computers per student, and expenditures per student. Determine whether the
# three variables have expanatory power by an F-Test of the hypothesis that all three
# have zero coefficients and via the Bayesian information criterion (BIC). The latter can
# be computed from an F-statistic: The BIC rejects the restriction when the F-statistic
# exceeds the log of the sample size.

reg1 <- "testscr ~ str + comp_stu + expn_stu"
lm1 <- lm(data = df, formula(reg1))
summary(lm1)

# (2) Do the same thing with a regression that adds the demographic variables: Average
# income, subsidized meals, calWorks per cent, and English learners percent. Again
# check whether the three “policy variables have explanatory power using an F test
# and BIC. Here you may need to extract the covariance matrix of coefficients from the
# lm()output to construct the F or chi-squared statistic.

reg2 <- paste0(reg1, " ")


# (3) Repeat the previous estimations and tests in models that add county fixed effects. In
# R using lm(), this is accomplished by just adding “county” to the list of right-hand
# side variables. (county is a “factor” in the R dataframe, so R automatically converts
# it into the appropriate array of dummy variables when including it in a regression.)



# (4) The districts vary greatly in size. Average scores might have more sampling variation
# in small districts. Plot the squared residuals from the estimated model in 2 against
# the total enrollment variable. Estimate a linear regression of these squared residuals
# on 1/enrl_tot. Use the inverse of these predicted values as the weights argument in 
# lm() (or otherwise estimate the corresponding weighted regression estimates) in
# the question 2 regression.

# (5) For at least two of the above regression models, calculate standard errors clustered
# by county. This is done very easily with the vcovCL() function from the sandwich
# package — so easily that if you’re doing it this way you might want to see how much
# difference it makes in all of the above regressions.

# (6) Estimate a random effects model, with county effects. In R, use the lme() function
# from the nlme package to estimate the 7-variable regression, with random effects by
# county. You do this by giving lme the argument random = ˜1 | county. Also use
# the argument method="ML", so that the estimation is by maximum likelihood.

# (7) Compare the random effects 7-variable model to the fixed effects model. In R, you
# can do this by re-estimating the fixed-effect model with the gls()function from the
# nlme package, again being sure to use method="ML" argument. 
# The summary() function applied to either random effects or fixed effects models 
# computed this way deliver both log likelihood and BIC values, so the models can be compared both by a
# frequentist chi-squared test based on the log likelihood and via the BIC.

# (8) Finally, for the random effects model, use a regression of its squared residuals on
# 1/(total enrollment) to generate weights for a weighted random effects estimation;
# see if this improves likelihood and/or changes important estimates. 
# [Note: I think that in the nlme estimation functions the “weights” arguments are variance scales —
# the inverse of the weights used in lm(). So you would use a weights= ˜w argument
# to lme() if you used weights=1/w in lm()].

# (9) Be ready to discuss: Does the evidence favor an important effect from the 
# “controllable” variables? The sizes and signs of the estimated effects, not just the significance
# levels of tests, should inform your views on this.


