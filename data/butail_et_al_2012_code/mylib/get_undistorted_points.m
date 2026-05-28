function [uv_undistorted uv_cost]=get_undistorted_points(uv, cam)
%function undistorted_points=get_undistorted_points(uv, cam)
%
% does an array search to find the undistorted points in the real world
% image plane. 
% uv is the 2x1 vector in pixels of the distorted point
% cam is the camera structure
% uv_undistored is a 2 x 1 vector in mm (or whichever units camera
% calibration is in)
% uv_cost is the distance between distorted pixel and the one found after
% mapping (should be less than .5 for any credibility)

ikm=inv(cam.km);

xp_ref=uv(1);
yp_ref=uv(2);


% to look within this region
x_range=xp_ref-100:1:xp_ref+100;
y_range=yp_ref-100:1:yp_ref+100;

[xp_ref yp_ref]=meshgrid(x_range, y_range);


% in mm
xr=ikm(1,1)*xp_ref + ikm(1,2)*yp_ref + ikm(1,3);
yr=ikm(2,1)*xp_ref + ikm(2,2)*yp_ref + ikm(2,3);


% distort these points
r2=xr.^2+yr.^2;
xd=xr.*(1+cam.kc1*r2+cam.kc2*r2.^2);
yd=yr.*(1+cam.kc1*r2+cam.kc2*r2.^2);

% back to pixels again
xp_d=xd*cam.km(1,1) + yd*cam.km(1,2) + cam.km(1,3);
yp_d=xd*cam.km(2,1) + yd*cam.km(2,2) + cam.km(2,3);

[r c]=size(xp_d);

% find your distorted reference pixel in this array
dist=sum((repmat(uv, 1, r*c)-[xp_d(:)'; yp_d(:)']).^2,1);
[val idx]=min(dist);

uv_undistorted=[xr(idx), yr(idx)]';

uv_cost=sqrt(val);