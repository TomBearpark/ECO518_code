%% Preliminaries

clear
clc
addpath('functions')
Raw = readtable('fish.csv');
df  = Raw;
df.C = df.logp*0+1;


XVars = {'C', 'logp'};
yVars = {'logq'};
RVars = {'C', 'stormy', 'mixed'};

X = df{:, XVars};
y = df{:, yVars};
R = df{:, RVars};

Res2SLS = est2SLS(X, R, y);

Res2Step = estIV2Step(X,R,y);
