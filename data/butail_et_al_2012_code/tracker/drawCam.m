function drawCam(fov, fr, cTw, frame_yn)
%function drawCam(fov, fr, cTw, frame_yn)
%
% fov is [2x1] field of view matrix fov(1) is angle of view wide in
% degrees, fov(2) is angle of view vertical in degrees
% fr is range of field of view. 
% cTw is the camera transformation matrix in the world frame. cTw is either
% [ 4 x 4 ] or [3 x 3]
% frame_yn is 0 if you don't want the orthogonal frame on the camera and 1 if you
% want
%
% Example:
% drawCam([30,45], 100, cTw, 1);


gca; hold on;

% if size(cTw,1) ==3  
%     cTw1(1:2,1:2) = cTw;
%     cTw1(1:3,3)= 
% exit

% camera cylinder radius
cylr=fr/25;
[camx camy camz] = cylinder(cylr);

% camera cylinder height
camz=camz*fr/10;

aovw = fov(1)*pi/180; % angle of view wide
aovv = fov(2)*pi/180;

% camera base side
side_b=fr/15;

xbase = [ -side_b, side_b, side_b, -side_b] ;
ybase = [ -side_b, -side_b, side_b, side_b];
zbase = [ 0 , 0 , 0, 0];

% creating field of view 
my=tan(aovv/2)*fr; mx=tan(aovw/2)*fr;
xfov=[0 mx mx; 0 -mx mx; 0 -mx -mx;0 mx -mx]';
yfov=[0 my -my;0 my my; 0 my -my; 0 -my -my]';
zfov=[0 fr fr; 0 fr fr; 0 fr fr; 0 fr fr]';

xyz=inv(cTw)*eye(4);
or=inv(cTw)*[0 0 0 1]';
if (frame_yn), plotFrame(or, xyz, fr/2, 3, eye(3),0); end
[camxg, camyg, camzg]=trpa2b(camx, camy, camz, inv(cTw));
surf(camxg, camyg, camzg,'FaceColor', 'k', ...
            'FaceAlpha', 1, 'EdgeColor', 'none'); hold on;
[xfovg yfovg zfovg]=trpa2b(xfov, yfov, zfov, inv(cTw));
patch(xfovg, yfovg, zfovg, [.5 .5 .5], 'EdgeColor', ...
            [.5 .5 .5], 'FaceAlpha',0.1);

% uncomment to create base for the camera        
% [xbaseg, ybaseg, zbaseg]=trpa2b(xbase, ybase, zbase, inv(cTw));
% patch(xbaseg, ybaseg, zbaseg, [.5 .5 .5], 'EdgeColor', ...
%             [.5 .5 .5], 'FaceAlpha',0.5);
