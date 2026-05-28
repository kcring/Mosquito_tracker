function mycube(l,w,h,c)
%function mycube(a,b,c,center)
% use patch to create a cube/cuboid with length(l), width(w), height(h),
% and center (c).
%
% c= [3x1]

c_mat=c*ones(1,4);
lwh_mat=[l;w;h]*ones(1,4)/2;

oper1=[  1 1 -1 -1;
        -1 1 1 -1;
        -1 -1 -1 -1];

oper2=[  1 1 -1 -1;
        -1 1 1 -1;
        1 1 1 1];
    

% front face
face(:,:,1)=(c_mat+lwh_mat.*oper1)';

% back
face(:,:,2)=(c_mat+lwh_mat.*oper2)';

% bottom
face(:,:,3)=(c_mat+lwh_mat.*circshift(oper1,2))';

% top
face(:,:,4)=(c_mat+lwh_mat.*circshift(oper2,2))';

% left
face(:,:,5)=(c_mat+lwh_mat.*circshift(oper1,1))';

% right
face(:,:,6)=(c_mat+lwh_mat.*circshift(oper2,1))';


gca; 
patch(squeeze(face(:,1,:)), squeeze(face(:,2,:)), squeeze(face(:,3,:)),  [.5 .5 .5], 'EdgeColor', ...
            [.5 .5 .5], 'FaceAlpha',0.1);
hold on;
