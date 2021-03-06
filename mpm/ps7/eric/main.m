%% Preliminaries
clear
clc

addpath('functions/')

figDir = 'figures/';
figSize = [6,3]; 
nGrid = 100;
scatteryLim = [0,1];
pVec = 1:10;


%% Specifications  
Res = struct;

yLab = 'sharefoodexp';
xLab = 'logincome';

% 1ii: Kernel regression
jSpec            = 1;
Res(jSpec).hGrid = 0:.01:1;
Res(jSpec).p     = 0;
Res(jSpec).yLab   = yLab;
Res(jSpec).xLab   = xLab; 
Res(jSpec).lab    = 'kernel';

% 1iii: Local linear regression
jSpec            = 2;
Res(jSpec).hGrid = 1:.01:2;
Res(jSpec).p     = 1;
Res(jSpec).yLab  = yLab;
Res(jSpec).xLab  = xLab; 
Res(jSpec).lab    = 'local linear';

% 1iv: Local linear regression

for p = pVec
    jSpec            = 2+p;
    Res(jSpec).hGrid = inf;
    Res(jSpec).p     = p;    
    Res(jSpec).yLab  = yLab;
    Res(jSpec).xLab  = xLab; 
    Res(jSpec).lab    = ['poly, p=' num2str(p)];
end

nSpec = length(Res);



%%

Raw              = readtable('engel.csv');
Raw.logincome    = log(Raw.income);
Raw.sharefoodexp = Raw.foodexp ./ Raw.income;



%% 1i: kernel density of log income
x         = Raw.logincome;
xGridKern = linspace(5, 9, nGrid);
F         = kerndensity(Raw.logincome, xGridKern);

% Plot density
f = figure;
resizeFig(figSize)
plot(xGridKern, F)
saveas(gcf, [figDir 'p1i_density.png'])



%% Run and plot specifications

for jSpec = 1:nSpec
    hGrid = Res(jSpec).hGrid;
    p     = Res(jSpec).p;
    y     = Raw{:, Res(jSpec).yLab};
    X     = Raw{:, Res(jSpec).xLab};
    lab   = Res(jSpec).lab;
    
    
    xGrid = linspace(min(X), max(X), 100);
    
    [mhat, CV, hOpt] = localpolyCV(X, y, xGrid, p, hGrid);
    Res(jSpec).mhat  = mhat;
    Res(jSpec).CV    = CV;
    Res(jSpec).hOpt  = hOpt;
    
    % Plot CV
    figure
    resizeFig(figSize)
    plot(hGrid, CV)
    title(['h^*_{CV}=' num2str(hOpt)])
    saveas(gcf, [figDir lab '_CV.png'])
    
    % Plot scatterplot and optimal regression curve
    figure
    resizeFig(figSize)
    plot(xGrid, mhat)
    hold on
    scatter(X, y)
    ylim(scatteryLim)
    saveas(gcf, [figDir lab '_scatter.png'])
    
end


%% Plot CV for polynomial

polyCV    = [Res(end-length(pVec)+1:end).CV];
[~, jOpt] = min(polyCV);

figure()
scatter(pVec, polyCV)
resizeFig(figSize)
xlabel('Order')
ylabel('CV')
box on
grid on
title(['p^*=' num2str(pVec(jOpt))])
saveas(gcf, [figDir 'poly_CV.png'])
ylim([0,.01])
saveas(gcf, [figDir 'poly_CV_tight.png'])




