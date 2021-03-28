function F = kerndensity(X, XI, h)
% Input
%   - X:  data
%   - XI: grid
% Output
%  - F: kernel density estimator
%
% Argument for bandwidth h is optional. If nothing is provided, default is
% the normal reference rule.

%%
X  = X(:);
XI = XI(:);

sigma = std(X);
N     = length(X);
NGrid = length(XI);


if ~exist('h')
    h = 1.059 * sigma / N^(1/5);
end



%%

XIMat = repmat(XI, 1, N);
XMat  = repmat(X', NGrid, 1); 



F = sum(normpdf((XMat - XIMat) /h), 2) / (N*h);

end