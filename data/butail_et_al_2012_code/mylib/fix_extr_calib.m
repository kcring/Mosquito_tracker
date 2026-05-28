function opt_trm= fix_extr_calib
% function auto_calibration
%
%

% Logic
% This function computes extrinsic calibration of multiple cameras given
% corresponding points between camera frames
% 1) Mark at least 30 corresponding points (np) between random camera frames chosen
% from the actual dataset (allow zooming in)
%
% 2) Use the initial camera calibration function to initialize the vector of
% unknown quantities. The unknown parameters are [R t], and [r_1, r_2,
% r_3]*np. 
%
% 3) Each pair of corresponding points gives us 3 new unknowns in the
% form of 3D position of that point and 6 repeating unknowns in the form of
% external parameters. 
%
% 4) The cost function to be minimized is sum((u_i-hat_u_i)^2 +
% (v_i-hat_v_i)^2) where hat_u_i, and hat_v_i are obtained by using the
% relation of perspective projection and zhang's distortion model

addpath ../external/TOOLBOX_calib/
addpath ../calib


floc='/home/sachit/Research/data/fish/2010/Aug/05/4fish03/';
get_cam_calib=@get_fishcam_calib_aug0510;

uvfile=[floc, 'output/data/', 'extrfix_uvpts.mat'];

% camids=['cam1'; 'cam2'; 'cam3'];
camids=['cam1'; 'cam2'];


tbf_camnum=2;

% Processing
nc=size(camids,1);

% Get list of all image files
for cam_ii=1:nc
    disp(sprintf('Listing files... ls %s%s*.*',  floc, camids(cam_ii,:)));
    c(cam_ii).flist = dir([floc, 'frames/', camids(cam_ii,:), '*.*']);
    if(~size(c(cam_ii).flist,1))
        error('Image files not found');
    end
    c(cam_ii).ni=size(c(cam_ii).flist,1);
end

% At least 30 randomly chosen points must be selected from randomly chosen
% images. After that more points 
total_points=12; % multiples of 6

% number of points clicked per image
npi=2;

% r0 for fmincon
rpts0=[100 100 -200];
rpts_range=[500 500 500];
cr_range=[20,20,20];
th_range=[pi/3 pi/3 pi/3];

clicked_colors=colormap(jet(npi));
ni=c(1).ni;


if(exist(uvfile, 'file'))
    load(uvfile);
    uv_jj=size(uv_clicked,1);
else
    uv_jj=0;
end

ni_choose=ceil(total_points/npi);

img_ids=unique(floor(rand(ni_choose,1)*(ni-1)))

img_jj=0;



figure(1); gcf; clf;
for img_ii=img_ids'
    img_jj=img_jj+1;
    proc=[];
    figure(1); gcf; 
    for cam_ii=1:nc
        subplot(1,nc, cam_ii); gca; cla;
        c(cam_ii).flist(img_ii).name
        img=imread([floc, 'frames/', c(cam_ii).flist(img_ii).name]);
        imshow(img); hold on;
    end
    
    proc=input('Image OK? Proceed? (zoom in if needed) ([]=yes, other=no):');
    if (isempty(proc))
        figure(1); gcf;
        for uv_ii=1:npi
            uv_jj=uv_jj+1;            
            for cam_ii=1:nc
                subplot(1,nc,cam_ii); gca;
                [uu vv]=ginput(1);
                plot(uu,vv, '+', 'Color', clicked_colors(uv_ii,:));
                uv_clicked(uv_jj, 2*cam_ii-1)=uu;
                uv_clicked(uv_jj, 2*cam_ii)=vv;
            end
        end
    end
end


% Minimization
n_points=size(uv_clicked,1);
n_eq=2*nc*n_points;
n_unknowns=6+3*n_points;

lsvar=n_eq-n_unknowns;
fprintf('# equations - # unknowns = %d\n', lsvar);

if(lsvar > 1)
    
    cam_tbf=get_cam_calib(tbf_camnum);

    disp('cTw0=');
    cam_tbf.trm
    th0=rodrigues(cam_tbf.trm(1:3,1:3));
    cr0=cam_tbf.trm(1:3,4);
    
    % **** Change this accordingly
    r0=repmat(rpts0', n_points,1);

    X0=[th0; cr0; r0];

    th_min=th0-th_range';
    cr_min=cr0-cr_range';
    rmin=repmat(-rpts_range', n_points,1);

    Xmin=[th_min; cr_min; rmin];

    th_max=th0+th_range';
    cr_max=cr0+cr_range';
    rmax=repmat(rpts_range', n_points,1);

    Xmax=[th_max; cr_max; rmax];

    for c_ii=1:nc
        camstrs(c_ii)=get_cam_calib(c_ii);
    end
    
    % constrained minimization
    options = optimset('Display','off','TolFun',1e-8);
    [opt_val, fval, residual]=fmincon(@(X) optfun(X, uv_clicked, camstrs, tbf_camnum), ...
                                            X0, ... 
                                            [],[],[],[],...
                                            Xmin,Xmax,...
                                            [], options);
    disp('Optimization done....')

    opt_rot_matrix=rodrigues(opt_val(1:3));

    opt_trm=[opt_rot_matrix, opt_val(4:6)
            0 0 0 1];
    
    X0(1:6)
    opt_val(1:6)
    
    residual

    fval

else
    fprintf('Need more points..\n');
end
save(uvfile, 'uv_clicked');
                                    
% optimization function                        
function fval=optfun(X, uv_clicked, camstrs, tbf_camnum)

fval=0;

n_points=size(uv_clicked,1);
nc=size(uv_clicked,2)/2;

tbf_cRw=rodrigues(X(1:3));
tbf_cTw=[tbf_cRw, X(4:6)
         0 0 0 1];
     
camstrs(tbf_camnum).trm=tbf_cTw;

for np_ii=1:n_points
    for c_ii = 1:nc
        pval=w2cam_nd(X(6+3*np_ii-2:6+3*np_ii), camstrs(c_ii));
        fval=norm(pval-uv_clicked(np_ii,2*c_ii-1:2*c_ii)')+fval;
    end
end
