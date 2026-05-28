function pts_ind=check_points_in_ellipse(ec, ea, eb, theta, pts)
%function pts_ind=check_points_in_ellipse(ec, ea, eb, theta, pts)
%
% function to check whether points pts exist within the ellipse boundary.
% This function can be used in ellipse fitting -- for example using least
% squares or for silhouette creation by snapping the pixels from within the
% bounding ellipse to the edge point.
%
% ec is a [2 x 1] center of ellipse 
% ea, eb are the major and minor axis of the ellipse
% theta is the orientation in radians
% pts is a [2 x n] matrix of n 2D points that must be checked
% image_yn is an indicator =1 for image coordinates and 0 for regular 2D (x
% to the right and y towards top)
% pts_ind is a [? x n] vector of indices for the points that lie within the
% ellipse

n=size(pts,2);

% transform them to the ellipse's frame
rotm=rota(theta, 'z');
rotm(1:2,3)=[ec(1), ec(2)]'; % we get the new frame in old frame

t_pts=inv(rotm)*[pts(1,:);pts(2,:); ones(1,n)];

tx=t_pts(1,:);
ty=t_pts(2,:);

tval=tx.^2/ea^2+ty.^2/eb^2-1;
pts_ind=find(tval<=0);
