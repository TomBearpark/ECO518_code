# Settings ----

library(sandwich)
library(lmtest)
library(tidyverse)
rm(list = ls(all.names = TRUE))
setwd('/Users/ericqian/Princeton/Courses/ECO 518/homework/hw6')
source('routines.R')


# DGP settings
alpha0  = 1
beta0   = 0
N       = 50
df      = 5


# Simulation settings
nMC   = 1000
nBoot = 500


Res     = sapply(1:nMC, function(x) runSim(x, nBoot, N, alpha0, beta0, df))

Table1 = matrix(nrow = nMC, ncol = 5)
Table2 = matrix(nrow = nMC, ncol = 5)
for(j in 1:5){
  Table1[,j] = sapply(Res['ex1',], '[[', j)
  Table2[,j] = sapply(Res['ex2',], '[[', j)
}



TableFinal = data.frame(parti = colMeans(Table1),
                        partii = colMeans(Table2))

row.names(TableFinal) = c('a', 'b', 'c', 'd', 'e')

write.csv(TableFinal, file = 'tableFinal.csv')

# Functions ---

