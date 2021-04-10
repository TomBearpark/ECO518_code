function out = makeTableProbit(Res, XVars)
% Part i
n               = Res.n;
out             = array2table([Res.beta sqrt(diag(Res.VMLE)/n) sqrt(diag(Res.VQMLE)/n)],...
    'VariableNames', {'Coefficient', 'SE_MLE', 'SE_QMLE'});
out.Properties.RowNames = XVars;
end