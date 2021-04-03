%% Preliminaries

clear
clc
addpath('functions')

Raw = readtable('fish.csv');
df  = Raw;

df.C = df.logp*0+1;


est2SLS(X, R, y)