function pc = get_plucker_coordinates(pts, cam)
% pts is [2 x n] in pixels
% cam is structure
% pc is [6 x n]

% plucker lines
wTc=inv(cam.trm);
iK=inv(cam.km);
n=tra2b([pts; ones(1, size(pts,2))], [wTc(1:3,1:3)*iK, [0 0 0]'; 0 0 0 1]);
n= normv(n);
m= cross(wTc(1:3,4)*ones(1,size(n,2)), n);
pc = [n; m];    