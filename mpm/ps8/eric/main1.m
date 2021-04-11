%% Preliminaries0+.05
clear
clc
addpath('functions')

Raw = readtable('mroz.csv');
df  = Raw;

% Add variables
df.C = df.age * 0 + 1;

%% 

covariates = {'C', 'kidslt6', 'age', 'educ', 'nwifeinc'};
y          = df.part;
X          = df{:, covariates};
Res        = estProbit(X ,y);

makeTableProbit(Res, covariates)


%% Test inversion

A = [strcmp(covariates, 'kidslt6') ;
    strcmp(covariates, 'educ')];

[Grid1, Grid2] = meshgrid(-1.25:.01:-.5, .05:.01:.25);
pVal           = nan(size(Grid1));
nGrid          = length(Grid1(:));

for j = 1:nGrid
    Res_j = wald(Res, A, [Grid1(j) Grid2(j)]');
    pVal(j) = Res_j.pW;
end

close all
contourf(Grid1, Grid2, pVal, 100, 'LineColor', 'flat')
hold on
view(2)
contour(Grid1, Grid2, pVal, [.05 .05], 'LineColor', 'r')
xlabel('kdslt6')
ylabel('educ')


 