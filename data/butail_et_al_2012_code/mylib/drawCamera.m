function drawCamera(fov, cam, cylr, range, frameyn)
%function drawCamera(fov, cam, cylr, range, frameyn)
%
% fov is [2x1] field of view matrix fov(1) is angle of view wide in
% degrees, fov(2) is angle of view vertical in degrees
% cam is the cam structure
% cylr is cylinder radius that represents the camera
% range is the depth of range shown with the camera in the same units as
% cyl radius
% frameyn is an indicator (0=don't show frame at camera center, 1= show frame at
% camera center)

gca; hold on;

[camx camy camz] = cylinder(1);

camx=camx*cylr; camy=camy*cylr; camz=camz*cylr*20;

aovw = fov(1)*pi/180; % angle of view wide
aovv = fov(2)*pi/180;

my=tan(aovv/2)*range; mx=tan(aovw/2)*range;
xfov=[0 mx mx; 0 -mx mx; 0 -mx -mx;0 mx -mx]';
yfov=[0 my -my;0 my my; 0 my -my; 0 -my -my]';
zfov=[0 range range; 0 range range; 0 range range; 0 range range]';

xyz=inv(cam.trm)*eye(4);
or=inv(cam.trm)*[0 0 0 1]';
if (frameyn), plotFrame(or, xyz, range/2); end
[camxg, camyg, camzg]=trpa2b(camx, camy, camz, inv(cam.trm));
surf(camxg, camyg, camzg,'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none'); hold on;
[xfovg yfovg zfovg]=trpa2b(xfov, yfov, zfov, inv(cam.trm));
patch(xfovg, yfovg, zfovg, [.5 .5 .5], 'EdgeColor', [.5 .5 .5], 'FaceAlpha',0.1);