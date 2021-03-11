#  A litle bit of code to go with solution to pset 6

library(dplyr)
library(xtable)

mat1 <- matrix(
  c(0.1, 0, 0.3, 0, 
    0, 0.9, 0, -0.6, 
    .6, 0, -0.3, 0, 
    0, .1, 0, .1), byrow = TRUE, nrow = 4
)
yt1 <- matrix(c(1,2,3,4))

mat1_tilde <- matrix(
  c(0.1, 0.3, 0, 0,
    .6,-0.3, 0, 0,
    0, 0, .9, -0.6, 
    0, 0, 0.1, .1), byrow = TRUE, nrow = 4
)
yt1_tilde <- matrix(c(1,3,2,4))

# Testwe get the same matrices
determinant(mat)
determinant(mat2)
# Check stationarity
eigen(mat)
mat1 %*% mat1

# Print coefficients in latex
mat1
mat1_tilde %*% mat1_tilde  %>% 
  xtable(align = rep("", 5), digits = 3) %>% 
  print(tabular.environment="bmatrix", digits = 3,
        hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
