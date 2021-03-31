function [mhat, CV, hOpt] = localpolyCV(X, y, X0, p, hGrid)

nh    = length(hGrid); 
CV    = nan(nh, 1);


for jh = 1:nh
    h         = hGrid(jh);
    
    % Print to console
    clc
    disp(['h=' num2str(h) '...'])
        
    [mHat, W] = localpoly(X, y, X, p, h);    
    CV(jh)    = mean(((y - mHat) ./ (1-diag(W)) ).^2);    
end

[~, jOpt] = min(CV);
hOpt      = hGrid(jOpt);

mhat = localpoly(X,y,X0,p, hOpt);
