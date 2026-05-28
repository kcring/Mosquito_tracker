function g = getrandSE3()

%{
    basis elements for se(3)
%}

E=zeros(4,4,6);

% rotation about x axis
E(:,:,1)=[  0 0 0 0;
            0 0 -1 0;
            0 1 0 0;
            0 0 0 0];

% rotation about y        
E(:,:,2)=[  0 0 1 0;
            0 0 0 0;
            -1 0 0 0;
            0 0 0 0];

% rotaiton about z
E(:,:,3)=[0 -1 0 0;
    1 0 0 0;
    0 0 0 0;
    0 0 0 0];

% along x
E(:,:,4)=[0 0 0 1;
    0 0 0 0;
    0 0 0 0;
    0 0 0 0];

% along y
E(:,:,5)=[0 0 0 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
% along z
E(:,:,6)=[0 0 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0];

stdb=[1 1 1 1 1 1]';
w(1)=-pi +pi*rand;
w(2)=-pi +pi*rand;
w(3)=-pi +pi*rand;
w(4:6)=rand(3,1);

X=zeros(4);
for i =1:6
    % ref: Murray et al. A mathematical introduction to Robotic
    % manipulation: pg 27, eq. 2.9
    X = X+stdb(i)*w(i)*E(:,:,i);
end

g= expm(X);

