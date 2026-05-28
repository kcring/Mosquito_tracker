function sd=stereo_matching(u1, u2, cam1, cam2, disp_limits)
%function sd=stereo_matching(u1, u2, cam1, cam2, disp_limits)
%
% u1 is [2x1] pixel position in camera 1 
% u2 is [2xn] pixel position in camera 2
% cam1 and cam2 are the camera structures
% disp_limits is a [2x1] vector of disparity limits [min max] (u1 is
% assumed to be the left camera i.e. u1>=u2)
%
% sd [nx1] is the distance between lines projecting out of the pixel points
% Good for a GNN type matching...
% 

n=size(u2,2);

% distorted positiion coordinates
% TBD: undistort (can be taken care of by passing undistorted coordinates)
p1=inv(cam1.km)*[u1;1];
p1=inv(cam1.trm)*[p1; 1]; p1=p1(1:3);
z1=inv(cam1.trm)*[0 0 0 1]'; z1=z1(1:3);
v1=p1-z1;
v1=v1/norm(v1);

u2(3,:)=1;
p2=inv(cam2.km)*u2;
p2(4,:)=1;
p2=inv(cam2.trm)*p2; p2=p2(1:3,:);
z2=inv(cam2.trm)*[0 0 0 1]'; z2=z2(1:3);
v2=p2-repmat(z2, 1, n);
v2=v2./repmat(sqrt(sum(v2.^2)),3,1);


sd=100*ones(n,1); % big number should be rejected outside
ind_valid=find(u1(1,:)-u2(1,:)>disp_limits(1) & u1(1,:)-u2(1,:)<disp_limits(2) & ...
                abs(u1(2,:)-u2(2,:))<20);
np_valid=size(ind_valid,2);
sd(ind_valid)=abs(dot(cross(repmat(v1,1,np_valid),v2(:,ind_valid)), repmat(p1, 1,np_valid)-p2(:,ind_valid)));