function [Xh_ P_]=initTarget(u, Xi)

Xh_=ones(Xi.nX,1);

Xh_(1:2,1)=u;
P_=diag([1 1 50 50]);
% Xh_(Xi.id,1)=id;
