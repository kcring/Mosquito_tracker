function dist=linePoint2dDist(pt, l)
% function dist=linePoint2dDist(pt, l)

% dist=abs(-l.m*pt(1,:)+pt(2,:)-l.c)/sqrt(1+l.m^2);
npt=size(pt,2);
dist=abs(dot(pt-l.r*ones(1,npt), [-l.u(2); l.u(1)]*ones(1,npt)));
