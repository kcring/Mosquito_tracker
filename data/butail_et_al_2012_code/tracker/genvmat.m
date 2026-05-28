function vmat=genvmat(cost, thresh)
% function vmat=genvmat(cost, thresh)
%

vmat=cost*0;
vmat(cost<thresh)=1;