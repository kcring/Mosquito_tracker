function vp=w2cam_rdot(tp,  cam, te, Xi)

ctp=[   tra2b(tp(Xi.ri,:), cam.trm)
        tra2b(tp(Xi.rdi,:), [cam.trm(1:3,1:3), zeros(3,1); 0 0 0 1])];

cnum1=ctp(Xi.rdi(1),:).*ctp(Xi.ri(3),:)-ctp(Xi.ri(1),:).*ctp(Xi.rdi(3),:);
cnum2=ctp(Xi.rdi(2),:).*ctp(Xi.ri(3),:)-ctp(Xi.ri(2),:).*ctp(Xi.rdi(3),:);
cden=ctp(Xi.ri(3),:).^2-ctp(Xi.rdi(3),:).^2*te^2/4;
        
% NOTE: multiplying by the time of exposure converts all this into the
% length terms instead of velocity terms
vn1=cnum1./cden*te;
vn2=cnum2./cden*te;

% ignoring distortion
% Note the 0 in the end, because we are projecting lines not points **
vp(1,:) = cam.km(1,1)*vn1 + cam.km(1,2)*vn2 + cam.km(1,3)*0;
vp(2,:) = cam.km(2,1)*vn1 + cam.km(2,2)*vn2 + cam.km(2,3)*0;