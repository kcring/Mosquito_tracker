
calibDate=input('Date when this calibration was done(mmmddyyyy):', 's');
extrinsicImageLocation='../calibration/';

fn_name=[extrinsicImageLocation, 'get_cam_calib_', calibDate, '.m'];

load('../calibration/camera_calibration.mat');

ff=fopen(fn_name,'w');
fprintf(ff,'\nfunction cam = get_cam_calib_%s(camid)', calibDate);
fprintf(ff,'\n%% ============================ ');
fprintf(ff,'\n%% Date of calibration:%s', calibDate);
fprintf(ff,'\n%% ============================ ');

fprintf(ff, '\n\nswitch camid');

for c_ii=1:size(cams,2)
    fprintf(ff,'\ncase %d', c_ii);
    fprintf(ff,'\ncam.id=\''%s\'';', cams(c_ii).id);
    fprintf(ff,'\ncam.km=[%f \t %f \t %f;\n%f \t %f \t %f;\n%f \t %f \t %f];',...
                cams(c_ii).fc(1), cams(c_ii).alpha_c*cams(c_ii).fc(1), cams(c_ii).cc(1), ...
                 0, cams(c_ii).fc(2), cams(c_ii).cc(2), ...
                 0, 0, 1);
             
    fprintf(ff,'\ncam.kc1=%f;', cams(c_ii).kc(1));
    fprintf(ff,'\ncam.kc2=%f;', cams(c_ii).kc(2));
    
    fprintf(ff,'\ncam.trm=[%f \t %f \t %f \t %f;', cams(c_ii).cTw(1,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f;', cams(c_ii).cTw(2,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f;', cams(c_ii).cTw(3,:));    
    fprintf(ff,'\n%f \t %f \t %f \t %f];', cams(c_ii).cTw(4,:));        

%     fprintf(ff,'\ncam.color=[%.3f %.3f %.3f];', cam_colors(c_ii,1), cam_colors(c_ii,2), cam_colors(c_ii,3));
end
fprintf(ff, '\nend');
fclose(ff);