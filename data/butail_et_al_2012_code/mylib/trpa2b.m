function [bx by bz] = trpa2b(ax, ay, az, bTa)
%function [bx by bz] = trpa2b(ax, ay, az, bTa)
%
% ax, ay, az can be any dimensions. 
% The output bx, by, bz are of the same dimensions. 
% bTa is a [4x4] transformation matrix
%
% See also TRA2B


bx = bTa(1,1).*ax+bTa(1,2).*ay+bTa(1,3).*az+bTa(1,4).*1;
by = bTa(2,1).*ax+bTa(2,2).*ay+bTa(2,3).*az+bTa(2,4).*1;
bz = bTa(3,1).*ax+bTa(3,2).*ay+bTa(3,3).*az+bTa(3,4).*1;
