%% Preliminaries
clear
clc
addpath('functions')

Raw = readtable('mroz.csv');
df  = Raw;

% Add variables
df.C = df.age * 0 + 1;

%% 



covariates = {'C', 'kidslt6', 'age', 'educ', 'nwifeinc'};
y = df.part;
X = df{:, covariates};
betaHat = estProbit(X ,y);