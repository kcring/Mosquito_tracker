function idx = resample(wts)
%function idx = resample(wts)
% wts:  weights (N x 1)
% idx:  sorted index values
%
% -SB, 2008

N=length(wts);
csw = cumsum(wts);
% start at bottom
i=1;
% draw a starting point
u1=rand*(1/N);
u=zeros(1,N);
for j = 1:N
    u(j) = u1+(1/N)*(j-1);
    while u(j) > csw(i)
        i=i+1;
    end
    idx(j)=i;
end