function fg=apply_threshold(fg, binary_t, area_t)
%function fg=apply_threshold(fg, binary_t, area_t)

fg=im2bw(fg,binary_t);
% fg=imfill(fg,'holes');

% fg=bwareaopen(fg,area_t, ones(3));
fg=bwareaopen(fg,area_t, [0 1 0; 1 1 1; 0 1 0]);
