function [eX eY] = myellipse(xc, yc, a, b, phi)
%function [eX eY] = myellipse(xc, yc, a, b, phi)
%
% xc, yc as the center
% a, b semi-major and minor axis
% phi the angle (in radians) semi-major axis makes with x axis. 

elem=50;
th=linspace(0,2*pi, elem);

x = sqrt(a^2)*cos(th);
y = sqrt(b^2)*sin(th);

trm=[cos(phi) -sin(phi) xc; 
        sin(phi) cos(phi) yc;
        0 0 1];

eX =  x*trm(1,1) + y*trm(1,2) + trm(1,3)*ones(1,elem);
eY =  x*trm(2,1) + y*trm(2,2) + trm(2,3)*ones(1,elem);