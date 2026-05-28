clc
clear all
close all
dr = '/media/puneetjain/Puneet_Back/dumpss/2011/';
imageNames = dir(fullfile(dr,'*.jpg'));
imageNames = {imageNames.name}';

outputVideo = VideoWriter(fullfile('/media/puneetjain/Puneet_Back/dumpss/','out_192011_2.avi'));
outputVideo.FrameRate = 30;
open(outputVideo)

for ii = 1:length(imageNames)
   img = imread(fullfile(dr,imageNames{ii}));
   writeVideo(outputVideo,img)
end
close(outputVideo)
