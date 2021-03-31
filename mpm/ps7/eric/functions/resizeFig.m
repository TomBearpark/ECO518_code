function resizeFig(dim)
% dim: [length, width]
f               = gcf;
f.Units         = 'inches';
f.Position(3:4) = dim;
end

