function val=normal2(x,mu,covar)
% x    : d x n; n data points with dim d
% mu   : d x 1
% covar: d x d covariance matrix
% Copyright 2009, Ali Bahramisharif
% This code is free to change, use and re-distribute.

[d, n]=size(x);
x = x-mu*ones(1,n);
x = x';

val = exp(-0.5*sum((x/(covar)).*x, 2))/ sqrt((2*pi)^d*abs((1e-10) +det(covar)));
