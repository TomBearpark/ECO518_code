zStar = 1.96;
aGrid = [-5:.01:5];
pi    = normcdf(aGrid-zStar) + normcdf(-aGrid - zStar);

plot(aGrid, pi)