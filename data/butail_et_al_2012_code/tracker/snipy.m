function c=snipy(cell, col, varargin)
%function r=snip1(cell, row, ri)
%
% snip1 gets you a section of large matrix made of cells
% row is the nth row of cells
% ri is the subindex within the cell

cc=cell(2)*(col-1)+1:cell(2)*col;

if nargin==3
    c=cc(varargin{1});
else
    c=cc;
end