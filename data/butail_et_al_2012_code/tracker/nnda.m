function [cost as Bjt] = nnda(costmat)
%function [cost as Bjt] = nnda(costmat)
%
% costmat is a [j x i] where j= # of measurements
% and i = # of targets + 1 (for clutter)
%
% cost is the total cost of matching

% sb Oct, 1, 2009
[nm nt]=size(costmat);

% Minimization is along the columns (for each target). In a
% scenario where the no. of measurements is 1, "[],1" is added so that it
% doesn't assume a vector
[cost as]=min(costmat, [], 1);
Bjt=zeros(nm,nt+1);
Bjt(:,1)=1;

for ii=1:nt
    if(as(ii))
        Bjt(as(ii),ii+1)=1;
        Bjt(as(ii),1)=0;
    end
end
cost=sum(cost);