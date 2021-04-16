%% Preliminaries
clear
clc
addpath('functions/')
figPath = 'figures/';


% Counterfactual settings for part 5
nBoot    = 1000;


Raw = readtable('heating.csv');
df  = Raw;
n = size(df, 1);


%% Organize data for input

% Choices crosswalk
choicesCW = table(string({'heatPump', 'gasCentral', 'electricCentral', ...
                'gasRoom', 'electricRoom'}'), [0 1 2 3 4]', ...
            'VariableNames', {'name', 'number'});

% choices Y (n x nChoices)
nChoices = height(choicesCW);
Y        = nan(n, nChoices);
for jChoice = 1:nChoices
    Y(:, jChoice) = df.choice == choicesCW.number(jChoice);    
end


% Covariates X (n x k x nChoices)
covariates = {'Ic', 'Oc'};  
X(:, 1, :) = df{:, append(choicesCW.name, 'Ic')};
X(:, 2, :) = df{:, append(choicesCW.name, 'Oc')};


Res      = estCondLogit(X, Y);
outLogit = table(Res.theta, Res.SE, 'VariableNames', {'thetaHat', 'SE'},...
    'RowNames', {'alpha1', 'alpha2', 'alpha3', 'alpha4', 'beta1', 'beta2'});

writetable(outLogit, [figPath, 'p2iii.xlsx']);

%% Part iv: Testing all constants jointly zero

ResNull = estCondLogit(X,Y, 'false');
LR      = 2*(Res.loglik - ResNull.loglik);
pLR     = 1 - chi2cdf(LR, nChoices-1); 

disp('LR statistic')
disp(LR)
disp('p value')
disp(pLR)

%% Part v: Market share

cfVar           = {'gasCentral', 'gasRoom'};
X1              = Res.X;
X2              = X1;

% Increase IC gas (index 6) cost for cfVar by 10%
X2(:, 6, [2,4]) = X2(:, 6, ismember(choicesCW.name, cfVar)) * 1.1;
delta           = CF(Res, X1, X2);

deltaBoot = nan(nBoot, nChoices);

rng(1)

% Bootstrap standard errors
for jBoot = 1:nBoot
    clc
    disp(['Iteration ' num2str(jBoot) ' of ' num2str(nBoot) '...'])
    keep  = randsample(1:n, n, 'true');
    X_j   = X(keep,:,:);
    Y_j   = Y(keep, :);
    Res_j = estCondLogit(X_j, Y_j, 'true', Res.theta);
    
    X1_j = Res_j.X;
    X2_j = X1_j;
    X2_j(:, 6, [2,4]) = X2_j(:, 6, ismember(choicesCW.name, cfVar)) * 1.1;
    
    deltaBoot(jBoot, :) = CF(Res_j, X1_j, X2_j);    
    
end

SEcf = std(deltaBoot);


outCF = table(delta(:), SEcf(:), 'VariableNames', {'coeff', 'SE'},...
    'RowNames', choicesCW.name);
writetable(outCF, [figPath, 'p2v.xlsx']);




