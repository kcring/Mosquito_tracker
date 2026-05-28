function r=snipx(cell, row, varargin)
%function r=snip1(cell, row, ri)
%
% snip1 gets you a section of large matrix made of cells
% row is the nth row of cells
% ri is the subindex within the cell

rc=cell(1)*(row-1)+1:cell(1)*row;

if nargin==3
    r=rc(varargin{1});
else
    r=rc;
end