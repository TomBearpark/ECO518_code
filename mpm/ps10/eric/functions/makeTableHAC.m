function out = makeTableHAC(Res)


Res.coeff

coeff = Res.coeff;
se    = Res.se;
pVal  = 2*(1-normcdf(abs(Res.coeff ./ Res.se)));

out = array2table([coeff, se, pVal], 'VariableNames', {'coeff', 'se', 'pVal'});
out.Properties.RowNames = ['Constant', Res.vars(1:end-1)];

end