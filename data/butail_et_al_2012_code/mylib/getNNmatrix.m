function NN = getNNmatrix(zhat, z, rgate)
% function NN = getNNmatrix(zhat, z, rgate)
%
% zhat is [z_k x nt ] matrix with columns equal to the number of targets
% z is a [z_k x nm ] measurement matrix with columns equal to the number of targets
% rgate is a metric in measurement space to validate measurements
% NN is a [nm x nt] matrix

% sb Oct, 1, 2009

[jk, nt]=size(zhat);
[jk, nm]=size(z);

% vm matrix with rows as measurements and columns as targets
NN=zeros(nm,nt); % nm x nt
invm=zeros(nm,nt); % matrix with innovations
        
% innovations are simply the absolute distance between ii-th
% measurement and jj-th target
for ii=1:nm
    for jj=1:nt
        invm(ii,jj)=norm(zhat(:,jj)-z(:,ii));
    end
end

% set the vm matrix entry to 1 for the smallest distance for that
% target
for jj=1:nm
    [tval tidx]=min(invm(jj,:));
    if(tval<rgate) % threshold or gate
        NN(jj,tidx)=1;
    end
end

% some targets will have multiple assignments. pick the min. out of those
for tt=1:nt
    [na jk]=size(find(NN(:,tt)>0));
    if(na>1)
        idxs=find(NN(:,tt)>0);
        [val pick]=min(invm(idxs,tt));
        for aa=1:na
            if(idxs(aa)~=pick)
                NN(idxs(aa),tt)=0;
            end
        end
    end
end 
