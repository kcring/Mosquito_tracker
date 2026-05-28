function cost=gencost(Z,zh, S)

cost=zeros(size(zh,2), size(Z,2));

% compute innovation 
for i=1:nz % measurements
    for j=1:nt % estimated measurements corresponding to targets
        inov=Z(:,i)-zh(:,j);
        cost(i,j)=inov'/S(:,:,j)*inov;
    end
end