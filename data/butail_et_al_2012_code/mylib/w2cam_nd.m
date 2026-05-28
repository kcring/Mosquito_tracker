function pix = w2cam_nd(rW, cam_str)
%function pix = w2cam(rW, cam_str)
%
% rW:       3 x wn matrix with world coordinates
% cam_str:  camera structure with cam.trm as the transformation matrix
% 
% pix:      2 x wn image coordinates
% ref: http://www.vision.caltech.edu/bouguetj/calib_doc/htmls/parameters.html


% get in camera frame
cr= tra2b(rW, cam_str.trm);

% normalized coordinates
xn = cr(1,:)./cr(3,:);
yn = cr(2,:)./cr(3,:);

% distorted coordinates
% rn2 =  xn.^2+yn.^2;
% xd = xn.*(1+cam_str.kc1*rn2 +cam_str.kc2*rn2.^2);
% yd = yn.*(1+cam_str.kc1*rn2 +cam_str.kc2*rn2.^2);
xd=xn; yd=yn;

% pixel coordinates
rp(1,:) = cam_str.km(1,1)*xd + cam_str.km(1,2)*yd + cam_str.km(1,3)*1;
rp(2,:) = cam_str.km(2,1)*xd + cam_str.km(2,2)*yd + cam_str.km(2,3)*1;

% return
pix=[rp(1,:);
    rp(2,:)];
