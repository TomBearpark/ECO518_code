function out = estCondLogit(Xin, Y, constant, theta0)
% estCondLogit()
% Input:
% This model contains a choice-dependent constant term.
%  - X: Covariates (n x nk x nOpts) for number of regressors nk and number
%       of options nOpts.
%  - Y: Dependent variable. 


arguments
   Xin(:,:,:) double  
   Y(:,:) double
   constant char = 'true'
   theta0 double = []
end


%% Setup
[n, nOpts] = size(Y);
k          = size(Xin, 2); 


% Setup constants
if strcmp(constant, 'true')
    nC         = nOpts-1;  % Number of constant terms
elseif strcmp(constant, 'false')
    nC = 0;
else
    error('Enter valid constant option...')
end

% If no initialization is provided
if isempty(theta0)
    theta0             = zeros(k+nC,1);
end


% Add columns corresponding to choice-dependent constant term. The first
% choice is normalized to zero by default. The constant terms indexed
% before the regressors
X = cat(2, zeros(n,nC, nOpts), Xin);
for jK = 1:nC
    X(:, jK, jK+1) = 1;
end


%% Maximum likelihood

fun                = @(theta) -loglikCondProbit(theta, X, Y);
options            = optimoptions('fmincon','Display', 'off');
[theta, negloglik] = fmincon(fun, theta0, [], [],  [], [], [], [], [], options);


%% Get standard errors

Omega = 0;  % Fisher information matrix

% Get information matrix
for i = 1:n
    X_i = squeeze(X(i,:,:));
    Y_i = squeeze(Y(i,:,:));
    
    p_i = exp(theta' * X_i)/sum(exp(theta' * X_i));    
    
    % alpha_0 is normalized to zero
    if strcmp(constant, 'true')
        alpha_i = Y_i(2:end) - p_i(2:end);
    else
        alpha_i = [];
    end
    beta_i  = (Y_i - p_i)*X_i(nC+1:end, :)';
    Delf    = [alpha_i beta_i]';
    Omega   = Omega+Delf*Delf';    
end

Omega = Omega/n;
V     = inv(Omega);  % Variance (under correct specification)
SE    = diag(sqrt(V/n));




%% Setup output

out          = struct;
out.theta    = theta;
out.loglik   = -negloglik;
out.V        = V;
out.SE       = SE;
out.constant = constant;
out.X        = X;
out.Y        = Y;
out.Xin      = Xin;
end