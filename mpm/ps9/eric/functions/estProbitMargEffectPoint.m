function mfx = estProbitMargEffectPoint(Res, XBar, jVar, varType, DeltaX)
% See estProbitMargEffect()


%% Preliminaries
beta  = Res.beta;  % Coefficients
X     = Res.X; 
n     = Res.n;

if ~strcmp(XBar, 'mean')
    XBar  = XBar(:);
end
%% Get marginal effect

if strcmp(varType, 'cts')  % cts marginal effect
    
    if strcmp(XBar, 'mean')
        
        mfx = mean(normpdf(X * beta) .* beta(jVar));
        
    else
        mfx = normpdf(XBar' * beta) * beta(jVar);
    end
    
elseif strcmp(varType, 'discrete')
    
    
    if strcmp(XBar, 'mean')
        [X1, X2]    = deal(X);
        X1(:, jVar) = DeltaX(1);
        X2(:, jVar) = DeltaX(2);
        
        
        mfx = mean(normcdf(X2 * beta)- normcdf(X1 * beta));
        
    else
        [XBar1, XBar2] = deal(XBar);
        XBar1(jVar)    = DeltaX(1);
        XBar2(jVar)    = DeltaX(2);        
        
        mfx = normcdf(XBar2' * beta) - normcdf(XBar1' * beta) ;
    end
    
    
    
else
    error('Invalid type.')
end


end
