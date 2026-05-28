function plotFrame(p, R, s, varargin)
%plotFrame(p, R, s, varargin)
%
% p is the position of the origin
% R is the rotation matrix
% s is the scale
%
% varargin{1} is lw, linewidth
% varargin{2} is color of each axis, row of [3x3] matrix 
% varargin{3} is whether arrow (1) should be used or quiver (0)

gca; hold on;

if(nargin>3)
    lw=varargin{1};
    color=varargin{2};
else
    lw=1.5;
    color=[ 1 0 0
            0 1 0
            0 0 0];
end

if nargin==6
    if varargin{3}
        arrow3(p(1:3)*ones(1,3), p(1:3)*ones(1,3)+R(1:3,1:3)*20, 'd-1', 1, 2);
%         arrow3([0,0,0; 0,0,0; 0,0,0], [200,0,0;0,200,0;0,0,500], 'd-1', 1, 2);
    end
else
    % x 
    quiver3(p(1), p(2), p(3), R(1,1), R(2,1), R(3,1), s, 'Color', color(1,:), 'LineWidth', lw);
    % y
    quiver3(p(1), p(2), p(3), R(1,2), R(2,2), R(3,2), s, 'Color', color(2,:), 'LineWidth', lw);
    % z
    quiver3(p(1), p(2), p(3), R(1,3), R(2,3), R(3,3), s, 'Color', color(3,:), 'LineWidth', lw);
    
end