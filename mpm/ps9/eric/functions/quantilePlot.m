function quantilePlot(Time, Quantiles, baseColor)
% Time:       Tx1  vector corresponding to x-margin
% Quantiles:  Txk: The k columns must be in increasing order along each
%             row. For example, Quantiles = [p5, p25, p50, p75 p90] is
%             valid. The center will be plotted as a solid line. Outer
%             quantiles are plotted with decreasing transparency.
% baseColor:  Gives the base band color.


% Turn time into a column vector
Time = Time(:);

if nargin < 3
    baseColor = [44,127,184]./255;  % Blue
    
end

if mod(size(Quantiles, 2), 2) ==0
    error('Enter an odd number of quantiles...')
end

nQ     = (size(Quantiles, 2)- 1)/2;
Center = Quantiles(:, nQ+1);
plot(Time, Center, 'LineWidth', 2, 'Color', baseColor*.9)
hold on

for jQ = 1:nQ
    bottom = Quantiles(:, nQ+1 - jQ);
    top    = Quantiles(:, nQ+1 + jQ);
    
    idx = ~isnan(bottom + top);
    fill([Time(idx); flip(Time(idx))], ...
        [bottom(idx); flip(top(idx))], baseColor,...
        'LineStyle','none',...
        'FaceAlpha', 1-jQ/(nQ+1));
end


end