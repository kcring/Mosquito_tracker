function [Xbe Ybe Zbe]= bentEllipsoid(a,b,c,K)
%function [Xbe Ybe Zbe]= bentEllipsoid(a,b,c,K)
%
% a,b,c are the lengths of semi-major,-medium, and minor axes
% K is the curvature coefficient

% Parameterize the ellipse coordinates
[Th Ph]=meshgrid(-pi:.1:pi, -pi:.1:pi);

% careful here when making the ellipse Note that a y-xbX^2 in the ellipse
% equation becomes a + xb*X^2 here
Xbe=a*sin(Th).*sin(Ph); 
Ybe=b*sin(Th).*cos(Ph) + K*Xbe.^2 ; 
Zbe=c*cos(Th);