function loglik = loglikCondProbit(theta, X, Y)
% loglikCondProbit()  Log likelihood function. 
theta = theta(:);

[n, ~] = size(Y);

loglik = 0;

for i = 1:n 
   X_i    = squeeze(X(i, :, :)); 
   Y_i    = Y(i, :);
   loglik = loglik + theta'*X_i*Y_i' - log(sum(exp(theta'*X_i)));
end

end