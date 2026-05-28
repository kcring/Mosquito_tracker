function R = completeOrthFrame(x)
%function R = completeOrthFrame(x)
%
% NOTE: make sure that the inertial frame is such that [0 0 1] => up
% x [3 X N] matrix of heading vector
% R [3 X 3 X N] rotation matrices

N=size(x,2);

vert=repmat([0,0,1]', 1,N);
y=cross(vert, x); 

% 2 norm
y_2=sqrt(y(1,:).^2+y(2,:).^2+y(3,:).^2);
y_2=repmat(y_2,3,1);
y=y./y_2;

z=cross(x,y); 
z_2=sqrt(z(1,:).^2+z(2,:).^2+z(3,:).^2);
z_2=repmat(z_2,3,1);
z =z./z_2; % z completes the orthonormal frame

for ii=1:N
    R(:,:,ii)=[x(:,ii), y(:,ii), z(:,ii)];
end
