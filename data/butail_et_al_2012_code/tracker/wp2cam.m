function [uC vC] = wp2cam(r1W,r2W,r3W,cam_str)
%function [uC vC] = wp2cam(r1W,r2W,r3W,cam_str)
%
% r1W,r2W,r3W:  same size matrices [m x n]
% cam_str:      camera structure with cam.trm as the transformation matrix
%
% uC, vC:       [m x n] sized outputs for image plane
% [r1W r2W r3W] to to cam_str camera coordinates. cam_str has a 

% ref: http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html


% get in camera frame
[r1C r2C r3C]= trpa2b(r1W,r2W,r3W, cam_str.trm);

% normalized coordinates
xn = r1C./r3C;
yn = r2C./r3C;

% distorted coordinates
rn2 =  xn.^2+yn.^2;
xd = xn.*(1+cam_str.kc1*rn2 +cam_str.kc2*rn2.^2);
yd = yn.*(1+cam_str.kc1*rn2 +cam_str.kc2*rn2.^2);

% pixel coordinates
uC = cam_str.km(1,1)*xd + cam_str.km(1,2)*yd + cam_str.km(1,3)*1;
vC = cam_str.km(2,1)*xd + cam_str.km(2,2)*yd + cam_str.km(2,3)*1;
