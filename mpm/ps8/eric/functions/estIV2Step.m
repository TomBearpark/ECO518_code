function Res2Step = estIV2Step(Z,X,y)


[n, kx] = size(X);
[~, kz] = size(Z);
Sxx     = X'*X/n;
Sxz     = X'*Z/n;
sxy     = X'*y/n;
G       = - X'*Z/n;      % partial g / partial delta = -z_i' x_i


Res2SLS = est2SLS(Z,X,y);
W       = inv(Res2SLS.Omega);  % Optimal weight matrix


delta = Sxz' * (W\Sxz) \ Sxz' *(W \ sxy);
eNew        = y - Z*delta;  % Residual
Omega       = X'* diag(eNew.^2) * X/n; % g_i*g_i' where g_i = (y_i-z_i'*delta)*x_i



V = (G'*W *G) \ G'* W * Omega * W*G * inv(G'*W*G);  % Heteroskedastic V

Res2Step         = struct;
Res2Step.Res2SLS = Res2SLS;
Res2Step.delta   = delta;
Res2Step.Omega   = Omega;
Res2Step.W       = W;
Res2Step.V       = V;


end