function [floc, frmloc]=getexpdirs

fprintf('Select the LOCATION of the EXPERIMENT DIR with subfolders data, calib...\n');
floc = uigetdir('.', 'LOCATION of the EXPERIMENT DIR with subfolders data, calib');
floc=[floc, '/'];
if exist([floc, filesep, 'frames'], 'dir')
    frmloc = [floc, filesep, 'frames', filesep];
    frames = dir([frmloc, '*.*']);
    if size(frames,1) < 10
        frmloc = uigetdir('.', 'LOCATION where the FRAMES are located');
    end
else
    frmloc = uigetdir('.', 'LOCATION where the FRAMES are located');
    frmloc = [frmloc, '/'];
end
