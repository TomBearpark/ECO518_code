%% Preliminaries

clear
clc

% A-R test parameters
alpha = .05;
bGrid = (-5:.01:1)';

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

AR_set = estAndersonRubin(X,R,y, bGrid, alpha);




%% Create output

Table2SLS  = makeTable(Res2SLS, XVars);
Table2Step = makeTable(Res2Step, XVars);


zCrit = abs(norminv(alpha/2));
[Table2Step{'logp', 'Coefficient'} - zCrit*Table2Step{'logp', 'SE'},...
    Table2Step{'logp', 'Coefficient'} + zCrit*Table2Step{'logp', 'SE'}]

Res2Step.pJ


% Check
%Res1 = fitlm(df, 'logp~stormy+mixed');
%df.logpHat = Res1.Fitted;
%fitlm(df, 'logq~logpHat')

