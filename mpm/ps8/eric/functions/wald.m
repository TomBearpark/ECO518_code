function ResWald = wald(Res, A, gamma)

beta = Res.beta;
V   = Res.VQMLE;
n   = Res.n;

r = size(A, 1);  % Number of restrictions

F  = n*(A*beta - gamma)'*((A*V*A') \ (A*beta - gamma));
pW = 1-chi2cdf(F, r);


% Output
ResWald     = struct;
ResWald.Res = Res;
ResWald.A       = A;
ResWald.gamma   = gamma;
ResWald.pW      = pW;



end