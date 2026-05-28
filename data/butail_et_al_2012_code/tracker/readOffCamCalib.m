function cams=readOffCamCalib(calibloc, nc)

% cams=strCam(nc);
addpath(calibloc);
calibfun=dir([calibloc, '/' '*calib*.m']);

if(size(calibfun,1))
    get_cam_calib=str2func(calibfun.name(1:end-2));
else
    error('Could not find calibration...');
end

for ii=1:nc
    cams(ii)=get_cam_calib(ii);
end
