function [optTrack, swarm_boundaries, camid_suffixes, ...
        Xi, dataloc, movieloc, calibloc, expname, ...
        image_id, nc, imageIds, get_cam_calib, cams]=initproc(floc, frmloc)


% Processing (do not change)
if ~strcmp(floc(end), '/')
    floc=[floc, '/'];
end

if ~exist(floc, 'dir')
    fprintf('[!] Could not find floc=%s\n',floc);
    return
end

Xi=strXi;

% assign frames location, data location etc.

if isempty(frmloc)
    frmloc=[floc, 'frames/'];
end

[optTrack, swarm_boundaries, camid_suffixes]=config;

calibloc=[floc, 'calib/'];
dataloc=[floc, 'output/data/'];
movieloc=[floc, 'output/movies/'];

[expname, image_id]=scan_expfile(floc);


nc=size(camid_suffixes,1);

imageIds=[[image_id, camid_suffixes(1)]; [image_id, camid_suffixes(2)]];

% Calibration information
addpath([floc '/calib']);
calibfun=dir([floc 'calib/' '*calib*.m']);
if(size(calibfun,1))
    get_cam_calib=str2func(calibfun.name(1:end-2));
else
    fprintf('[?] Could not find calibration for the cameras. Is this 2D tracking? ...\n');
end

cams=readOffCamCalib(calibloc, nc);

warning('off', 'Images:initSize:adjustingMag');
warning('off', 'optim:fmincon:SwitchingToMediumScale');

