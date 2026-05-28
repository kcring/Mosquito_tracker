function [fval r] = localize_lsq(r0, m, get_cam_calib, cams, rmin, rmax)
%function [fval r] = localize_lsq(r0, m, get_cam_calib, cams, rmin, rmax)
%
% r0 is 3 x 1 matrix
% m is 2 x nc where nc is the number of cams in cams
% get_cam_calib is the address to calibration function 
% cams is the list of cams for the measurements m
% rmin, rmax are each 3 x 1 matrices that can 

% set optimization options
options = optimset('Display','off','TolFun',1e-6);

[r, fval, exitflag, output]=fmincon(@(r) optfun(r, m, get_cam_calib, cams), ...
                r0, ... 
                [],[],[],[],...
                rmin,rmax,...
                [], options);
          
fprintf('[I] _fmincon_ algorithm=%s\n',output.algorithm);
fprintf('[I] function calls=%d\n', output.funcCount);

        
        
function f=optfun(r, m, get_cam_calib, cams)
f=0;

for c = cams' % hardcoding the number of measurements each time
    cam=get_cam_calib(c);
    pix=w2cam(r, cam);
    xp=pix(1); yp=pix(2);
    f=sqrt((xp-m(1,c))^2+(yp-m(2,c))^2)+f;
end
        