%% Preliminaries 
clear
clc
addpath('functions/')
outDir = 'output/';

% Read data
Raw = readtable('jtpa.csv');
df  = Raw;

%%

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


