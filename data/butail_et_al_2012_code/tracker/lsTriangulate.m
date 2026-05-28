function [r, err]= lsTriangulate(Uk, cams)
%function [r err]= lsTriangulate(Uk, cams)
%
% Uk is 2xnc matrix


% use top and side views to get pose
nc=size(Uk,2);
con=zeros(nc*2,4);
for zz=1:nc
    % least squares triangulation
    con(2*(zz-1)+1,:)= [cams(zz).km(1,1), cams(zz).km(1,2), ...
        cams(zz).km(1,3)-Uk(1,zz)]*cams(zz).trm(1:3,1:4);
    con(2*(zz-1)+2,:)= [cams(zz).km(2,1), cams(zz).km(2,2), ...
        cams(zz).km(2,3)-Uk(2,zz)]*cams(zz).trm(1:3,1:4);
end

hh_A= con;
b= -hh_A(:,4);
A= hh_A(:,1:3);


r=(A'*A)\A'*b;

err=zeros(2,nc);
for zz=1:nc
    reproj=cams(zz).km*cams(zz).trm(1:3,1:4)*[r;1];
    err(:,zz)=Uk(:,zz)-reproj(1:2)/reproj(3);
end
err=norm(err);

