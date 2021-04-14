clear
clc
addpath('functions/')
addpath('../../ps8/eric/functions/')  % PS8 function path
figPath = 'figures/';

Raw = readtable('mroz.csv');
df  = Raw;

% Add variables
df.C = df.age * 0 + 1;

%% 

covariates = {'C', 'kidslt6', 'age', 'educ', 'nwifeinc'};
y          = df.part;
X          = df{:, covariates};
Res        = estProbit(X ,y);


jVar = find(strcmp(covariates, 'educ'));
XBar = mean(X);
mfx1 = estProbitMargEffect(Res, XBar, jVar);
mfx2 = estProbitMargEffect(Res, 'mean', jVar);

jVar = find(strcmp(covariates, 'kidslt6'));
mfx3 = estProbitMargEffect(Res, 'mean', jVar, 'discrete', [0, 1]);

%% Organize results
labs      = {'ME at mean', 'APE educ', 'APE kidslt6'};
out       = table(labs', 'VariableNames', {'parameter'});
out.value = [mfx1.mfx, mfx2.mfx, mfx3.mfx]';
out.SE    = [mfx1.SE, mfx2.SE, mfx3.SE]';

writetable(out, [figPath, 'p1.xlsx'])
