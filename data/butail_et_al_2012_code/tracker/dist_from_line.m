function dist = dist_from_line(pl, pts)
% function dist = dist_from_line(pl, pt)
%
% to compute distance from line
% pl is the plucker coordinates of a line (n,m) [6x1]
% where n is [3 x 1] vector representing the line direction
% m is also [3 x 1] moment of any point p on that line computed by pxn
% pts is point(s) in 3D space [3 x n] 
%
% dist is the same dim as pt
%
% ref: Rosenhahn 2007, "Three-dimensional shape knowledge for joint image
% segmentation and pose tracking"

n=pl(1:3);
m=pl(4:6);

nhat_T= [0 n(3) -n(2); -n(3) 0 n(1); n(2) -n(1) 0];

pxn(1,:)=pts(1,:)*nhat_T(1,1) + pts(2,:)*nhat_T(1,2) + pts(3,:)*nhat_T(1,3);
pxn(2,:)=pts(1,:)*nhat_T(2,1) + pts(2,:)*nhat_T(2,2) + pts(3,:)*nhat_T(2,3);
pxn(3,:)=pts(1,:)*nhat_T(3,1) + pts(2,:)*nhat_T(3,2) + pts(3,:)*nhat_T(3,3);

dist=sum((pxn-m*ones(1,size(pts,2))).^2);
