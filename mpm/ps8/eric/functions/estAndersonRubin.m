function [interval, rej, p] = estAndersonRubin(X, R, y, bGrid, alpha)

[n, k] = size(X);

assert(all(X(:, 1) == 1))  % Make sure first column is a one vector (constant)
assert(all(R(:,1) == 1))

rej = nan(n, 1);
p   = nan(n,1);

nGrid = length(bGrid);

Restr = [0 1 0;
         0 0 1];


for jb = 1:nGrid
    b   = bGrid(jb, :);
    e = y - X(:, 2:end)*b;
    
    % Use HC1 SE, consistent with STATA
    [estCov, ~, coeff] = hac(R(:, 2:end), e, 'type', 'HC', 'weights',...
   'HC1','display','off');    

    % Rk: MATLAB uses different notation from Hayashi. r is a restriction
    % function where r(theta)=0 on parameters of the restricted model.
    % R is the Jacobian of the restriction function evaluated at the
    % unrestricted model parameter.
    [rej(jb), p(jb)] = waldtest(coeff(2:end), Restr, estCov, alpha);
    
end

interval = quantile(bGrid(~rej), [0,1]);



end