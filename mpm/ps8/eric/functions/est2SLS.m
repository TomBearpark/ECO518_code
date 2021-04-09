function Res = est2SLS(Z, X, y)
% est2SLS()  Estimate 2SLS. Use Hayashi notation.
% Input:
%  - Z endogenous variables
%  - X instruments
%  - y dependent variable

% Setup
[n, kx] = size(X);
[~, kz] = size(Z);
G       = - X'*Z/n;      % partial g / partial delta = -z_i' x_i
W       = inv(X'*X/n);        % Weight matrix under 2SLS
Sxx     = X'*X/n;
Sxz     = X'*Z/n;
sxy     = X'*y/n;

% Estimate and get variance covariance
delta = (Sxz' * (Sxx \ Sxz)) \ Sxz' * (Sxx \ sxy);
e     = y - Z*delta;  % Residual
Omega = X'* diag(e.^2) * X/n; % g_i*g_i' where g_i = (y_i-z_i'*delta)*x_i
V     = (G'*W *G) \ G'* W * Omega * W*G * inv(G'*W*G);  % Heteroskedastic V


% Organize output object
Res       = struct;
Res.delta = delta;
Res.V     = V;
Res.Omega = Omega;
Res.X     = X;
Res.y     = y;
Res.Z     = Z;
Res.n     = n;
end