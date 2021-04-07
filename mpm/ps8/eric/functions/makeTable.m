function out = makeTable(Res, XVars)
% Part i
n               = Res.n;
out             = array2table([Res.delta sqrt(diag(Res.V)/n)], 'VariableNames', {'Coefficient', 'SE'});
out.Properties.RowNames = XVars;
out.tStat       = out.Coefficient./out.SE;
end