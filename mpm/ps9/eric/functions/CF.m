function delta = CF(Res, X1, X2)
% CF()  Counterfactual for conditional probit. Reports average difference
%       between probabilities implied by X1 and X2
% Input: 
%  - Res: Results from conditional probit estimation
%  - X1: Data
%  - X2: Data

X     = Res.X;
n     = size(X,1);
nOpts = size(X,3);
theta = Res.theta;
p1 = nan(n, nOpts);  % Unconditional
p2 = nan(n, nOpts);  % Counterfactual

for i = 1:n
    X1_i = squeeze(X1(i,:,:));
    X2_i = squeeze(X2(i,:,:));    
    
    p1(i,:) = exp(theta' * X1_i)/sum(exp(theta' * X1_i));    
    p2(i,:) = exp(theta' * X2_i)/sum(exp(theta' * X2_i));

end

delta = mean(p2 - p1);

end