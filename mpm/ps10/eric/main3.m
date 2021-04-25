%% Preliminaries

clear
clc
addpath('functions/')
addpath('../../ps8/eric/functions/')
outDir = 'output/';


yob = 50:53;  % DOB of interest
covariateLabel = {'constant', 'veteran'};



%% (i) Get IV coefficients and standard errors

Raw = readtable('draft.csv');
df  = Raw;

Spec        = struct;

for j = 1:length(yob)
    
    % Run IV for specific year
    Spec(j).yob  = yob(j);            
    keep         = df.yob == yob(j);
    Spec(j).keep = keep;
    
    
    % Setup matrices
    Spec(j).X    = [ones(sum(keep), 1) df.draftelig(keep) ];  % Instrument
    Spec(j).y    = df.lwage(keep);                            % Dependent
    Spec(j).Z    = [ones(sum(keep), 1) df.veteran(keep)];     % Endogenous variables

    % Run IV
    Spec(j).Res      = est2SLS(Spec(j).Z, Spec(j).X, Spec(j).y);        
    Spec(j).outTable = makeTable(Spec(j).Res, covariateLabel);
    
end


%% Format output

% Organize coefficients
coeffTable = cell2mat(cellfun(@(x) x.Coefficient, {Spec.outTable}, 'UniformOutput', false));
coeffTable = array2table(coeffTable, 'RowNames', covariateLabel, ...
    'VariableNames', append("yr", string(yob)));
writetable(coeffTable, [outDir, 'p3_coeff.xlsx'])


% Organize standard errors
SETable    = cell2mat(cellfun(@(x) x.SE, {Spec.outTable}, 'UniformOutput', false));
SETable = array2table(SETable, 'RowNames', covariateLabel, ...
    'VariableNames', append("yr", string(yob)));
writetable(coeffTable, [outDir, 'p3_SE.xlsx'])


%% Get fraction of compliers, always-takers, and never-takers

for j = 1:length(yob)
    
    % Get cohort
    df_j = df(Spec(j).keep, :);
    
    Spec(j).alwaysTakers = ...
        sum((df_j.veteran & (1-df_j.draftelig))) / sum((1-df_j.draftelig));  % veteran and draft ineligible
    Spec(j).compliers = ...
        sum((df_j.veteran & df_j.draftelig)) / sum((df_j.draftelig))...  % veteran and draft ineligible
        -Spec(j).alwaysTakers;
    Spec(j).neverTakers = 1 - Spec(j).compliers-Spec(j).alwaysTakers;
end


%% Would expect draft to be stronger in the general population. Never takers are a large 
% portion. Those who never take may have information on earnings power in
% the counterfactual where they enter draft (do what they can to avoid). As
% a result

idTable = array2table([[Spec.alwaysTakers]; [Spec.compliers]; [Spec.neverTakers]], ...
        'VariableNames', append("yr", string(yob)), ...
        'RowNames', {'alwaysTakers', 'compliers', 'neverTakers'});
    
writetable(idTable, [outDir, 'p3_id.xlsx'])
