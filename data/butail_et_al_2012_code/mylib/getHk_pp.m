function Hk = getHk_pp(r, cam)
%function Hk = getHk_pp(r, cam)
%
% r is 3x1 vector of world coordinates
% cam is the camera structure

if numel(r)>3, r=r(1:3); end

% convert to camera coordinates and normalize
cr=cam.trm*[r;1];
un=cr(1)/cr(3);
vn=cr(2)/cr(3);
rn2=un^2+vn^2; rn4=rn2^2;

% partial camera with world coordinates
d_cr1_r1=cam.trm(1,1);
d_cr1_r2=cam.trm(1,2);
d_cr1_r3=cam.trm(1,3);

d_cr2_r1=cam.trm(2,1);
d_cr2_r2=cam.trm(2,2);
d_cr2_r3=cam.trm(2,3);

d_cr3_r1=cam.trm(3,1);
d_cr3_r2=cam.trm(3,2);
d_cr3_r3=cam.trm(3,3);

% parital normalized with camera coordinates
d_un_cr1=1/cr(3);
d_un_cr2=0;
d_un_cr3=-cr(1)/cr(3)^2;

d_vn_cr1=0;
d_vn_cr2=1/cr(3);
d_vn_cr3=-cr(2)/cr(3)^2;

% partial distorted with normalized coordinates
d_ud_un=1+cam.kc1*rn2 + cam.kc2*rn4 + un*(2*cam.kc1*un + 4*cam.kc2*un*rn2);
d_ud_vn=un*(2*cam.kc1*vn + 4*cam.kc2*vn*rn2);

d_vd_un=vn*(2*cam.kc1*un + 4*cam.kc2*un*rn2);
d_vd_vn=1+cam.kc1*rn2 + cam.kc2*rn4 + vn*(2*cam.kc1*vn + 4*cam.kc2*vn*rn2);


% partial pixel with distorted coordinates
d_up_ud=cam.km(1,1);
d_up_vd=cam.km(1,2);

d_vp_ud=cam.km(2,1);
d_vp_vd=cam.km(2,2);



% Now start going back one at a time

d_un_r1=d_un_cr1*d_cr1_r1 + d_un_cr2*d_cr2_r1 + d_un_cr3*d_cr3_r1;
d_un_r2=d_un_cr1*d_cr1_r2 + d_un_cr2*d_cr2_r2 + d_un_cr3*d_cr3_r2;
d_un_r3=d_un_cr1*d_cr1_r3 + d_un_cr2*d_cr2_r3 + d_un_cr3*d_cr3_r3;

d_vn_r1=d_vn_cr1*d_cr1_r1 + d_vn_cr2*d_cr2_r1 + d_vn_cr3*d_cr3_r1;
d_vn_r2=d_vn_cr1*d_cr1_r2 + d_vn_cr2*d_cr2_r2 + d_vn_cr3*d_cr3_r2;
d_vn_r3=d_vn_cr1*d_cr1_r3 + d_vn_cr2*d_cr2_r3 + d_vn_cr3*d_cr3_r3;


d_ud_r1=d_ud_un*d_un_r1 + d_ud_vn*d_vn_r1;
d_ud_r2=d_ud_un*d_un_r2 + d_ud_vn*d_vn_r2;
d_ud_r3=d_ud_un*d_un_r3 + d_ud_vn*d_vn_r3;


d_vd_r1=d_vd_un*d_un_r1 + d_vd_vn*d_vn_r1;
d_vd_r2=d_vd_un*d_un_r2 + d_vd_vn*d_vn_r2;
d_vd_r3=d_vd_un*d_un_r3 + d_vd_vn*d_vn_r3;


d_up_r1=d_up_ud*d_ud_r1 + d_up_vd*d_vd_r1;
d_up_r2=d_up_ud*d_ud_r2 + d_up_vd*d_vd_r2;
d_up_r3=d_up_ud*d_ud_r3 + d_up_vd*d_vd_r3;


d_vp_r1=d_vp_ud*d_ud_r1 + d_vp_vd*d_vd_r1;
d_vp_r2=d_vp_ud*d_ud_r2 + d_vp_vd*d_vd_r2;
d_vp_r3=d_vp_ud*d_ud_r3 + d_vp_vd*d_vd_r3;

Hk=[d_up_r1 d_up_r2 d_up_r3 0 0 0 
    d_vp_r1 d_vp_r2 d_vp_r3 0 0 0];