function rotm = rota(a, axis)
%function rotm = rota(a, axis)
% a:    angle in radians nx1
% axis: 'x','y', or 'z'
%
% rotm: rotation matrix 3x3x n

switch axis
    case 'x'
        rotm=[1 0 0;
                0 cos(a) -sin(a);
                0 sin(a) cos(a)];
    case 'y'
        rotm=[cos(a) 0 sin(a);
                0 1 0;
                -sin(a) 0 cos(a)];
    case 'z'
        rotm=[cos(a) -sin(a) 0;
                sin(a) cos(a) 0;
                0 0 1];
    otherwise
        error('huh!');
end