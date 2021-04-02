%% Preliminaries

clear
clc
addpath('functions/')
figDir = 'figures_problem2/';
figSize = [6,3]; 

nBoot = 250;

Raw = readtable('nls.csv');

df                          = Raw;
df.Properties.VariableNames = {'Y', 'X', 'R'};
df.R2                       = df.R.^2;
df.C                        = df.R .^0;
n = size(df, 1);

lm1 = fitlm(df, 'Y~X+R+R2');



%% Point estimates

% E(Y|R), NW 
%hGrid_Y_R = [.1:.1:.7, .75:.025:1.25, 1.3:.1:1.5];
hGrid_Y_R                  = .75:.025:1.25;
[Y_RHat, CV_Y_R, hOpt_Y_R] = localpolyCV(df.R, df.Y, df.R, 0, hGrid_Y_R);


% E(X|R), NW
hGrid_X_R = .7:.01:.8;
[X_RHat, CV_X_R, hOpt_X_R] = localpolyCV(df.R, df.X, df.R, 0, hGrid_X_R);

beta0  = regress( df.Y - Y_RHat, [df.X - X_RHat (df.X - X_RHat)*0+1]);  




%% Non-parametric bootstrap

rng(1)
ResBoot = array2table(nan(nBoot,2), 'VariableNames', {'linear', 'semipara'});


for jBoot = 1:nBoot
    clc
    disp([num2str(jBoot) ' of ' num2str(nBoot) '...'])
    dfBoot = datasample(df, n, 'Replace', true);

    % Estimate LM
    lm1                   = fitlm(dfBoot, 'Y~X+R+R2');
    ResBoot.linear(jBoot) = lm1.Coefficients{'X', 'Estimate'};
    
    % Estimate semi parametric model
    [Y_RHatBoot] = localpoly(dfBoot.R, dfBoot.Y, dfBoot.R, 0, hOpt_Y_R);
    [X_RHatBoot] = localpoly(dfBoot.R, dfBoot.X, dfBoot.R, 0, hOpt_X_R);
    ResBoot.semipara(jBoot) = ...
        regress( dfBoot.Y - Y_RHatBoot, dfBoot.X - X_RHatBoot);
end



%% Make plots

% Bandwidth for NW E(Y|R)
figure()
plot(hGrid_Y_R, CV_Y_R)
resizeFig(figSize)
title(['h^{*}_{CV} = ' num2str(hOpt_Y_R)])
box on; grid on;
saveas(gcf, [figDir, 'CV_ExpYR.png'])


% Bandwidth for NW E(X|R)
figure()
plot(hGrid_X_R, CV_X_R)
resizeFig(figSize)
title(['h^{*}_{CV} = ' num2str(hOpt_Y_R)])
box on; grid on;
saveas(gcf, [figDir, 'CV_ExpXR.png'])

% Estimate for beta0 for NW 
figure()
scatter(df.X - X_RHat, df.Y - Y_RHat)
resizeFig(figSize)
title(['\beta_0=' num2str(beta0,3)])
saveas(gcf, [figDir, 'CV_scatterStage2.png'])


figure()
resizeFig([6,4])
hold on
ksdensity(ResBoot.semipara)
ksdensity(ResBoot.linear)
box on
grid on
lgd = legend({['Semi-parametric (\beta=' num2str(beta0, 3) ...
    ', SE=' num2str(std(ResBoot.semipara),3) ')'],...
    ['Series (p=2) (\beta=' num2str(lm1.Coefficients{'X', 'Estimate'}, 3) ...
    ', SE=', num2str(std(ResBoot.linear), 3) ')' ]}, ...
    'Location', 'southoutside');
saveas(gcf, [figDir, 'bootCompare.png'])

