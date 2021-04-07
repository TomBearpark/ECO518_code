function Res = estIV2Step(Z,X,y)


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

Res         = struct;
Res.Res2SLS = Res2SLS;
Res.delta   = delta;
Res.Omega   = Omega;
Res.W       = W;
Res.V       = V;
Res.X       = X;
Res.y       = y;
Res.Z       = Z;
Res.n       = n;

end