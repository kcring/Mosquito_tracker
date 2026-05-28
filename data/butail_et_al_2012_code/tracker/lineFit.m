function p=lineFit(x,y)
%function p=quadCurveFit(x,y)
%
% fits a quadratic curve to x,y
%
% y=p(1) + p(2)*x + p(3)*x.^3

np=numel(x);

A=[ones(np,1), x];


b=y;

if numel(y) > 3
    p=(A'*A)\A'*b;
else
    p=[0 1];
end
