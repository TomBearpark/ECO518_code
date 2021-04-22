%% Preliminaries

clear
clc
addpath('functions/')
addpath('../../ps8/eric/functions/')

yob = 50:53;  % DOB of interest
covariateLabel = {'constant', 'veteran'};


%%

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

coeffTable = cell2mat(cellfun(@(x) x.Coefficient, {Spec.outTable}, 'UniformOutput', false));
coeffTable = array2table(coeffTable, 'RowNames', covariateLabel, ...
    'VariableNames', append("yr", string(yob)));

SETable    = cell2mat(cellfun(@(x) x.SE       , {Spec.outTable}, 'UniformOutput', false));
SETable = array2table(SETable, 'RowNames', covariateLabel, ...
    'VariableNames', append("yr", string(yob)));





