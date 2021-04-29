%% Preliminaries 
clear
clc
addpath('functions/')
outDir = 'output/';

% Read data
Raw = readtable('jtpa.csv');
df  = Raw;

%% Run regressions

ResNoControls = struct;
ResNoControls.vars = {'treatment', 'earnings'};
% No controls
[ResNoControls.EstCov, ResNoControls.se, ResNoControls.coeff] = ...
    hac(df(:, ResNoControls.vars), 'type', 'HC', ...
    'weights', 'HC1');   % for consistency with STATA

outNoControls = makeTableHAC(ResNoControls);

writetable(outNoControls, [outDir 'p2_noControls.xlsx'])

% Base category is 55-78
ResControls = struct;
controls = {'age2225', 'age2629', 'age3035', 'age3644', 'age4554'};
ResControls.vars = ['treatment', controls, 'earnings'];
[ResControls.EstCov, ResControls.se, ResControls.coeff] = ...
    hac(df(:, ResControls.vars), 'type', 'HC', ...
    'weights', 'HC1');

outControls = makeTableHAC(ResControls);
writetable(outControls, [outDir 'p2_controls.xlsx'])


%% Get sample standard deviations of earnings for treatment/control


sigmaTreatment = std(df.earnings(df.treatment == 1));
sigmaControl   = std(df.earnings(df.treatment == 0));
a              = 1000;
propTreat      = 2/3;
propControl    = 1-propTreat;
zStar         = 1.96;  % size 5% test

gamma = @(n)  a*sqrt(n)/sqrt(sigmaTreatment^2/propTreat + sigmaControl^2/propControl);
pi    = @(n) normcdf(gamma(n) - zStar) + normcdf(-zStar - gamma(n));
 
 
nOpt = fmincon(@(n) abs(pi(n) - 0.8), height(df));
ceil(nOpt)
