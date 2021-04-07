%% Preliminaries

clear
clc
addpath('functions')
Raw = readtable('fish.csv');
df   = Raw;
df.C = df.logp*0+1;


XVars = {'C', 'logp'};
yVars = {'logq'};
RVars = {'C', 'stormy', 'mixed'};

X = df{:, XVars};
y = df{:, yVars};
R = df{:, RVars};


%% Results

Res2SLS  = est2SLS(X, R, y);
Res2Step = estIV2Step(X,R,y);


%% Create output

makeTable(Res2SLS, XVars)
makeTable(Res2Step, XVars)

% Check
%Res1 = fitlm(df, 'logp~stormy+mixed');
%df.logpHat = Res1.Fitted;
%fitlm(df, 'logq~logpHat')

