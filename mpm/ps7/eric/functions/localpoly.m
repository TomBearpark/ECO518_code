function [mHat, W] = localpoly(X, y, X0, p, h)
% Input:
%  - X:  independent variable
%  - y:  dependent variable 
%  - X0: grid points
%  - p:  Degree of polynomial
%  - h:  bandwidth
% Output
%  - m: fitted value
%  - W: weight matrix


[n, k] = size(X);

if k > 1
   assert('Not written for k>1!') 
end

X0  = X0(:);
nx0 = size(X0, 1); 

% Structure matrices. Z = [1, X, X^2,...], z = [1, X0, X0^2,...] 
Z = ones(n,1);
z = ones(nx0, 1);

for jp = 1:p
    Z = [Z, X.^jp];
    z = [z, X0.^jp];
end

% Matrix versions of X and evaluation points X0 (n by nx0)
XX  = repmat(X  , 1, nx0);
XX0 = repmat(X0', n, 1); 
KK  = normpdf((XX - XX0)/h);
MM  = nan(nx0, p+1);  % M = z'*(sum K((Xj-x)/h)Z_jZ_j')^-1

% Precompute M for speed
for j0 = 1:nx0
    MM(j0, :) = ((Z'* diag(KK(:, j0))* Z) \ z(j0,:)')';
end

W    = MM*Z' .* KK';  % Weight matrix
mHat = W * y;       % Predicted value

% W(50,:)*y
% 
% % 
%  scatter(X, y)
%  hold on
%  plot(X0, mHat)
end