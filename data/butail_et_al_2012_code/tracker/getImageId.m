function [frmloc image_id expname]=getImageId
% function [frmloc image_id]=getImageId
% 
% works with images from streampix with format XYZ_HH.MM.SS.####

% get the image id automatically

[filename, pathname]=uigetfile('*.bmp; *.tif; *.png', 'Choose a file in the sequence');
frmloc=[pathname, '/'];

tmp_name=textscan(filename, '%s', 'Delimiter', '.');
tmp_name=tmp_name{:};
expname=sprintf('%s.', tmp_name{1:end-2});
expname=expname(1:end-1);

tmp_name=textscan(filename, '%s', 'Delimiter', '_');
tmp_name=tmp_name{:};
image_id=sprintf('%s_', tmp_name{1:end-2});