function Res = estProbit(X, y)
% estProbit()
% Input
%  - X: Independent variables (n x k)
%  - y: Dependent variables (n x 1)

[n, k] = size(X);
assert(length(y) == n)

% Negative log likelihood function
negloglik = @(beta) -sum(y .* log(normcdf(X*beta)) + (1-y) .* log(1 - normcdf(X*beta)));
beta0     = (X'*X) \ X'*y ;  % Initialize with OLS

options = optimoptions('fmincon','Display', 'off');
betaHat  = fmincon(negloglik, beta0, [], [],  [], [], [], [], [], options);

% Compute weight
z = 2*y-1;  % Transform y to {-1, 1}. Easier PDF to handle

w     = normpdf(z.* X*betaHat) .*z ./ normcdf(z.* X*betaHat);
H     = -X'*diag((X*betaHat + w) .* w) * X/n;
Omega = X' * diag(w.^2) * X/n;
invH  = inv(H);

% Variance-covariance matrices.
VMLE  = -invH;
VQMLE = invH*Omega * invH;

% Output objects
Res       = struct;
Res.beta  =  betaHat;
Res.VMLE  = VMLE;
Res.VQMLE = VQMLE;
Res.n     = n;
Res.X     = X;
Res.y     = y;
end