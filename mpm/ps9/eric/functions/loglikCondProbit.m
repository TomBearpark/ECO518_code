function loglik = loglikCondProbit(theta, X, Y)
% loglikCondProbit()  Log likelihood function. 
theta = theta(:);

[n, k, nOpts] = size(X);


% Naive loop: 4x slower than vectorized version
% loglik = 0;
% for i = 1:n 
%    X_i    = squeeze(X(i, :, :)); 
%    Y_i    = Y(i, :);
%    loglik = loglik + ...
%             theta'*X_i*Y_i' ...   % First term
%     - log(sum(exp(theta'*X_i)));  % Second term
% end

% Vectorized version. See above code reference
loglik =...
    theta'*sum(squeeze(... % First term
    sum(X .* permute(repmat(Y, 1, 1, k), [1 3 2]), 3)))' ...  
     -sum(...  % Second term: sum across n observations
    log(...
    sum(...  
    exp(...
    squeeze(sum(... theta'*X_i
    permute(repmat(theta, 1, nOpts, n), [3 1 2]) .* X, 2))), 2)));  % Reorganize



end