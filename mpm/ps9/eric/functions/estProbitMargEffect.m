function out = estProbitMargEffect(Res, XBar, jVar, varType, DeltaX, SEType, nBoot)
% Input
%  - Res:    Results file from probit estimation procedure
%  - XBar:   Marginal effect evaluated at XBar. Use 'mean' to specify an
%            average marginal effect.
%  - jVar:   Variable of interest for marginal effect (index).
%  - type:   Discrete or continuous. Continuous is the default.
%  - DeltaX: For discrete marginal effect. Specify the change. (x_1, x_2).
%            Use the notation of the lecture notes:
%                G(x_2^j beta_0^j +...) - G(x_1^j beta_0^j +...)
%%
n = Res.n;
X = Res.X;
y = Res.y;

if ~exist('varType', 'var')
   varType = 'cts';
end

if strcmp(varType, 'discrete')
    if ~exist('DeltaX', 'var')
        error('Specify change in independent variable.')
    end
else
    DeltaX = [];
end



%% Get marginal effect
mfx = estProbitMargEffectPoint(Res, XBar, jVar, varType, DeltaX);


%% Get standard errors
    if ~exist('SEType', 'var')
        SEType = 'boot';
    end

    if ~exist('nBoot', 'var')
        nBoot = 1000;
    end


if strcmp(SEType, 'boot')
    MfxBoot = nan(nBoot, 1);
    
    % Non-parametric bootstrap
    for jBoot = 1:nBoot
        
        clc
        disp(['Running ' num2str(jBoot) ' of ' num2str(nBoot) ' bootstrap iterations...'])

        
        keep           = randsample(1:n, n, true);
        Res_j          = estProbit(X(keep, :), y(keep));
        MfxBoot(jBoot) = estProbitMargEffectPoint(Res_j, XBar, jVar, varType, DeltaX);
    end
    
    SE = std(MfxBoot);
end

%% Output

out        = struct;
out.mfx    = mfx;
out.SE     = SE;
out.SEType = SEType;
out.Res    = Res;
out.XBar   = XBar;
out.jVar   = jVar;
out.type   = varType;
out.DeltaX = DeltaX;
out.nBoot  = nBoot;

end