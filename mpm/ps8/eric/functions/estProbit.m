function betaHat = estProbit(X, y)
% estProbit()
% Input
%  - X: Independent variables (n x k)
%  - y: Dependent variables (n x 1)

[n, k] = size(X);
assert(length(y) == n)

% Negative log likelihood function
negloglik = @(beta) -sum(y .* log(normcdf(X*beta)) + (1-y) .* log(1 - normcdf(X*beta)));
beta0     = (X'*X) \ X'*y ;  % Initialize with OLS
betaHat   = fmincon(negloglik, beta0);


end