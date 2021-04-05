function out = est2SLS(Z, X, y)
% est2SLS()  Estimate 2SLS. Use Hayashi notation.
% Input:
%  - Z endogenous variables
%  - X instruments
%  - y dependent variable

[n, kx] = size(X);
[~, kz] = size(Z);
G       = - X'*Z/n;      % partial g / partial delta = -z_i' x_i
W       = X'*X/n;        % Weight matrix under 2SLS
Sxx     = X'*X/n;
Sxz     = X'*Z/n;
sxy     = X'*y/n;

delta = (Sxz' * (Sxx \ Sxz)) \ Sxz' * (Sxx \ sxy);
e     = y - Z*delta;  % Residual

Omega = X'* diag(e.^2) * X/n; % g_i*g_i' where g_i = (y_i-z_i'*delta)*x_i

V = (G'*W *G) \ G'* W * Omega * W*G * inv(G'*W*G);  % Heteroskedastic V


% Organize output object
out       = struct;
out.delta = delta;
out.V     = V;
out.Omega = Omega;

end