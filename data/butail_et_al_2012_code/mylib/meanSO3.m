function mR = meanSO3(rots,wts)
% Get mean rotation as per Moakher
% Usage: meanSO3(rots, weights) where rots are 3x3 matrices with associated
% weights
% Returns: 3x3 mean rotation

[jk,jk, N] = size(rots);

bR=zeros(3);
for i =1:N
    bR = wts(i)*rots(:,:,i) + bR;
end

[U S Vt]=svd(bR');
V = Vt';

if (det(bR') > 0)
    mR=V*U';
else
    mR=V*diag([1,1,-1])*U';
end
