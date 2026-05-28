function F= get_F_for_stereo(get_cam_calib)
% function F = get_F_for_stereo(get_cam_calib)
%
% this function assumes that the camera origin is at camera1
% ref: page 244, Multiple view geometr

cam1=get_cam_calib(1);
cam2=get_cam_calib(2);

% Pl=cam1.km*[eye(3), zeros(3,1)];
% Pr=cam2.km*cam1.trm(1:3,1:4);

t=cam2.trm(1:3,4);

% epipole of the first camera
% e=cam1.km*cam2.trm(1:3,1:3)'*t;

% cross product operator
% ex=[0 -e(3) e(2); e(3) 0 -e(1); -e(2) e(1) 0];
tx=[0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0];

% F=cam2.km'\cam2.trm(1:3,1:3)*cam1.km'*ex;
E=tx*cam2.trm(1:3,1:3); % essential matrix
F=cam2.km'\E/cam1.km;

% sanity checks
% fprintf('||Fe||=%.2f\n', norm(F*e));
% fprintf('||F^Tep||=%.2f\n', norm(F'*cam2.km*t));
% fprintf('rank(F)=%d\n', rank(F));