% generate camera calib function
% this script will generate a camera calibration function for a day using
% the extrinsic calibration from the images specified. It will also show a
% mean and variance of the extrinsic parameters
%
% The following functions are called from the calibration toolbox
% extract_grid, compute_extrinsic

clear all;

addpath ../external/TOOLBOX_calib;

%% change these values

% location where intrinsic calibration .mat files are kept
intrlocs=[  '/home/sachit/Research/calib/fishcams/feb21_2011/cam1/';
            '/home/sachit/Research/calib/fishcams/feb21_2011/cam2/';
            '/home/sachit/Research/calib/fishcams/feb21_2011/cam3/'];
            
% corresponding images for extrinsic. the identifiers are used to list the
% image
camid_list=['20110221_cam1'; '20110221_cam2'; '20110221_cam3'];        
extrimg_loc='/home/sachit/Research/calib/fishcams/feb21_2011/extr/run1/';

    
% processing
nc=size(camid_list, 1);
cam_colors=colormap(lines(size(camid_list,1)));

calibdate=input('Date when this calibration was done(mmmddyy):', 's');

fn_name=strcat('/tmp/get_cam_calib_', calibdate, '.m');


for ii=1:nc
       
    matfile=strcat(intrlocs(ii,:), 'Calib_Results.mat');
    load(matfile, 'fc','cc','kc','alpha_c');
    
      
    % compute the extrinsic for all images available for this id
    flist=dir(strcat(extrimg_loc, camid_list(ii,:),'*.*'));
    ni=size(flist,1);
    
    for im_ii=1:ni
        disp(strcat('processing ', flist(im_ii).name, '....'));
        I=double(imread(strcat(extrimg_loc, flist(im_ii).name)));
        if size(I,3)>1,
            I = I(:,:,2);
        end;
        wintx = 5;
        winty = 5;
        
        % function calls from camera calibration toolbox
        [x_ext,X_ext,n_sq_x,n_sq_y,ind_orig,ind_x,ind_y] = extract_grid(I,wintx,winty,fc,cc,kc);
        [omc_ext,Tc_ext,Rc_ext,H_ext] = compute_extrinsic(x_ext,X_ext,fc,cc,kc,alpha_c);
        
%       for mirror reflections do one inversion
%       careful: the new camera frame will be left-handed
%         if (ii==1)
%             Rc_ext=Rc_ext*[ 0 1 0
%                             1 0 0
%                             0 0 1];
%                        
%         end

        [x_reproj] = project_points2(X_ext,omc_ext,Tc_ext,fc,cc,kc,alpha_c);
        err_reproj = x_ext - x_reproj;
        err_std2 = std(err_reproj')';
        
        
        fprintf(1,'\n\nExtrinsic parameters:\n\n');
        fprintf(1,'Translation vector: Tc_ext = [ %3.6f \t %3.6f \t %3.6f ]\n',Tc_ext);
        fprintf(1,'Rotation vector:   omc_ext = [ %3.6f \t %3.6f \t %3.6f ]\n',omc_ext);
        fprintf(1,'Rotation matrix:    Rc_ext = [ %3.6f \t %3.6f \t %3.6f\n',Rc_ext(1,:)');
        fprintf(1,'                               %3.6f \t %3.6f \t %3.6f\n',Rc_ext(2,:)');
        fprintf(1,'                               %3.6f \t %3.6f \t %3.6f ]\n',Rc_ext(3,:)');
        fprintf(1,'Pixel error:           err = [ %3.5f \t %3.5f ]\n\n',err_std2); 

        
        cam(ii).extr_calib(:,:,im_ii)=[Rc_ext, Tc_ext;
                    0 0 0 1];
    end
    
    
    % this part is to display the results
    Basis = [X_ext(:,[ind_orig ind_x ind_orig ind_y ind_orig ])];

    VX = Basis(:,2) - Basis(:,1);
    VY = Basis(:,4) - Basis(:,1);

    nX = norm(VX);
    nY = norm(VY);

    VZ = min(nX,nY) * cross(VX/nX,VY/nY);

    Basis = [Basis VZ];

    [x_basis] = project_points2(Basis,omc_ext,Tc_ext,fc,cc,kc,alpha_c);

    dxpos = (x_basis(:,2) + x_basis(:,1))/2;
    dypos = (x_basis(:,4) + x_basis(:,3))/2;
    dzpos = (x_basis(:,6) + x_basis(:,5))/2;



    figure(2);
    image(I);
    colormap(gray(256));
    hold on;
    plot(x_ext(1,:)+1,x_ext(2,:)+1,'r+');
    plot(x_reproj(1,:)+1,x_reproj(2,:)+1,'yo');
    h = text(x_ext(1,ind_orig)-25,x_ext(2,ind_orig)-25,'O');
    set(h,'Color','g','FontSize',14);
    h2 = text(dxpos(1)+1,dxpos(2)-30,'X');
    set(h2,'Color','g','FontSize',14);
    h3 = text(dypos(1)-30,dypos(2)+1,'Y');
    set(h3,'Color','g','FontSize',14);
    h4 = text(dzpos(1)-10,dzpos(2)-20,'Z');
    set(h4,'Color','g','FontSize',14);
    plot(x_basis(1,:)+1,x_basis(2,:)+1,'g-','linewidth',2);
    title('Image points (+) and reprojected grid points (o)');
    hold off;

    print('-dpng', sprintf('%s/extr_cam_%d.png', extrimg_loc, ii));


    cam(ii).fc=fc;
    cam(ii).cc=cc;
    cam(ii).kc=kc;
    cam(ii).alpha_c=alpha_c;
    
    clear fc cc kc alpha_c
end

% reset the world frame to the first camera frame
for c_ii=nc:-1:1
    for im_ii=1:ni
        cam(c_ii).extr_calib(:,:,im_ii)=cam(c_ii).extr_calib(:,:,im_ii)*inv(cam(1).extr_calib(:,:,im_ii));
    end
end

% find the average of all extrinsic calibrations
for c_ii=1:nc
    if ni>1
        std(cam(c_ii).extr_calib,0,3)      
        for im_ii=1:ni
            rvec(:,im_ii)=rodrigues(cam(c_ii).extr_calib(1:3,1:3,im_ii));
        end
        
        cam(c_ii).cTw=[rodrigues(mean(rvec,2)), mean(cam(c_ii).extr_calib(1:3,4,:),3);
                                    0       0       0       1];
    else
        cam(c_ii).cTw=cam(c_ii).extr_calib;
    end
end


ff=fopen(fn_name,'w');
fprintf(ff,'\nfunction cam = get_cam_calib_%s(camid)', calibdate);

fprintf(ff,'\n%% ============================ ');
fprintf(ff,'\n%% Date of calibration:%s', calibdate);
fprintf(ff,'\n%% ============================ ');

fprintf(ff, '\n\nswitch camid');

for c_ii=1:nc
    fprintf(ff,'\ncase %d', c_ii);
    fprintf(ff,'\ncam.id=\''%s\'';', camid_list(c_ii,:));
    fprintf(ff,'\ncam.km=[%f \t %f \t %f;\n%f \t %f \t %f;\n%f \t %f \t %f];',...
                cam(c_ii).fc(1), cam(c_ii).alpha_c*cam(c_ii).fc(1), cam(c_ii).cc(1), ...
                 0, cam(c_ii).fc(2), cam(c_ii).cc(2), ...
                 0, 0, 1);
             
    fprintf(ff,'\ncam.kc1=%f;', cam(c_ii).kc(1));
    fprintf(ff,'\ncam.kc2=%f;', cam(c_ii).kc(2));
    
    fprintf(ff,'\ncam.trm=[%f \t %f \t %f \t %f;', cam(c_ii).cTw(1,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f;', cam(c_ii).cTw(2,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f;', cam(c_ii).cTw(3,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f];', cam(c_ii).cTw(4,:));        

    fprintf(ff,'\ncam.color=[%.3f %.3f %.3f];', cam_colors(c_ii,1), cam_colors(c_ii,2), cam_colors(c_ii,3));
end
fprintf(ff, '\nend');
fclose(ff);
    
