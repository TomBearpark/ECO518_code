%% Preliminaries
clear
clc
addpath('functions/')

figPath = 'figures/';

%
nBoot      = 500;
alpha      = .05;
qGrid      = .05:.01:.95;
figDim     = [6,6];
covariates = {'Constant', 'educ', 'exper', 'exper2'};

%%

% Read/clean data, setup objects
Raw         = readtable('nls.csv');
df          = Raw;
df.exper2   = df.exper.^2;
df.Constant = df.educ*0+1;
n           = height(df);
nq          = length(qGrid);
k           = length(covariates);
Beta        = nan(k, nq);          % Full sample
BetaBoot    = nan(k, nq, nBoot);   % Bootstrap


% Full sample
X = df{:, covariates};
y = df.luwe;


% Non parametric bootstrap
for jBoot = 1:nBoot
    clc
    disp(['Running ' num2str(jBoot) ' of ' num2str(nBoot) ' bootstrap draws...'])
    
    for jq = 1:nq
        
        
        if jBoot == 1
            Beta(:, jq) = rq(X,y,qGrid(jq));
        end
        
        % Make bootstrap dataset
        keep = randsample(n, n, true);
        X_j  = X(keep, :);
        y_j = y(keep);
        
        % Run quantile regression
        BetaBoot(:, jq, jBoot) = rq(X_j, y_j, qGrid(jq));
        
    end
end


%% Clean/output figure


zStar         = -norminv(alpha/2);
BetaBootSE    = std(BetaBoot, 0, 3);
BetaCI        = nan(size(Beta, 1), size(Beta,2),3);
BetaCI(:,:,1) = Beta - BetaBootSE .* zStar;
BetaCI(:,:,2) = Beta;
BetaCI(:,:,3) = Beta + BetaBootSE .* zStar;


close all
for jVar = 1:k
    subplot(2,2,jVar)
    quantilePlot(qGrid, squeeze(BetaCI(jVar, :, :)))
    xlabel('quantile')
    box on
    grid on
    title(covariates{jVar})
end
resizeFig(figDim)
saveas(gcf, [figPath 'p3_quantiles.png'])
