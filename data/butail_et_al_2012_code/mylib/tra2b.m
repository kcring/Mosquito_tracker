function rB = tra2b(rA, bTa)
%function rB = tra2b(rA, bTa)
%
% Transforms point(s) rA from frame A to frame B through
% matrix transformation bTa. rA is [2 or 3]xn. bTa is 3x3 or 4x4. 
% The output rB is of the same dimension as rA.
%
% See also TRPA2B

% - SB Jan 15, 2009
% changed Jun 12, 2010 to include planar matrices

if(size(bTa,1)-1 ~= size(rA,1))
    error('[!] incompatible matrices');
end

if(size(bTa,1)==4)
    xB = bTa(1,1).*rA(1,:)+bTa(1,2).*rA(2,:)+bTa(1,3).*rA(3,:)+bTa(1,4).*1;
    yB = bTa(2,1).*rA(1,:)+bTa(2,2).*rA(2,:)+bTa(2,3).*rA(3,:)+bTa(2,4).*1;
    zB = bTa(3,1).*rA(1,:)+bTa(3,2).*rA(2,:)+bTa(3,3).*rA(3,:)+bTa(3,4).*1;

    rB=[xB
        yB
        zB];
elseif(size(bTa,1)==3)
    xB = bTa(1,1).*rA(1,:)+bTa(1,2).*rA(2,:)+bTa(1,3).*1;
    yB = bTa(2,1).*rA(1,:)+bTa(2,2).*rA(2,:)+bTa(2,3).*1;

    rB=[xB
        yB];
end
    