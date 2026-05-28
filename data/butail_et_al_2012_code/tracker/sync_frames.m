% set these numbers to where both videos are synchronous
lssync=1;
hssync=1;

% if verify=1 then it will create a montage, otherwise it will create a
% dataset for tracking
verify=1; 

% if the frames are not in sync at the first frame itself
if lssync*2 < hssync
    hssync=hssync-2*lssync;
    lssync=1;
end

pref='EXP_2016_03_25'; % date of experiment (keep exact same format)
suff='19.17.17'; % time of experiment (keep exact same format)

% low speed camera
lsdir='/Users/sachit/Sandbox/mq_marker_exp/30fps_deinterlaced/';
lsfiles=dir([lsdir, '/*.bmp']);
lsfiles=lsfiles(lssync:end); 

% high speed camera
hsdir='/Users/sachit/Sandbox/mq_marker_exp/documents-export-2016-04-13/';
hsfiles=dir([hsdir, '/*.bmp']);
hsfiles=hsfiles(hssync:end);

% where to dump the frames
mndir='/Users/sachit/Sandbox/mq_marker_exp/20160325_191717/frames/';

for ii=1:size(hsfiles, 1)

    imls=imread([lsdir, lsfiles(ceil(ii/2)).name]);
    imhs=imread([hsdir, hsfiles(ii).name]);
    
    if verify
        % if you want a montage
        img1=[imls(1:1024,:,:), imhs];
        imwrite(img1, [mndir, sprintf('montage_%.5d.bmp', ii)]);
    else
        % else just write frames
        imwrite(imls(1:1024,:,:), [mndir, sprintf('%s_R_%s%.5d.bmp', pref, suff, ii)]);
        imwrite(imhs, [mndir, sprintf('%s_L_%s%.5d.bmp', pref, suff,  ii)]);
    end
    
end
    