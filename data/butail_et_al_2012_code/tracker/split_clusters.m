function CL= split_clusters(k, CL, Xh_, P_, Zk, MHT, Kalman, Xi)

for cl=1:size(CL,2)
    cln=1;
    
    if k-MHT.nscanback >0 && ~isempty(CL(cl).scan(k-MHT.nscanback).hyp)
    hypmat=cat(1,CL(cl).scan(k-MHT.nscanback).hyp.assignment);
    utest=zeros(size(hypmat,2),1);
    for ii=1:size(hypmat,2)
        col=unique(hypmat(:,ii));
        col=col(col~=0);
        utest(ii)=numel(col);
       
    end
    tid_split=hypmat(1,utest==1);
    
    
    
    for tt=tid_split(1:end-1)
        split=1;
        for kk=k-MHT.nscanback+1:k-1
            
            hypmat=cat(1,CL(cl).scan(kk).hyp.assignment);
            utest=zeros(size(hypmat,2),1);
            for jj=1:size(hypmat,2)
                utest(jj)=numel(unique(hypmat(:,jj)));
            end 
            if ~ismember(tt, hypmat(1,utest==1))
                split=0;
            end
        end
        if split
            tid(cln)=tt;
            zid(cln)=find(hypmat(1,:)==tt);
            cln=cln+1;
        end
    end
    ncl=size(CL,2);
    for cc=1:cln-1
        for kk=k-MHT.nscanback:k-1

            for jj=1:size(CL(cl).scan(kk).hyp,2)
                ii=find(CL(cl).scan(kk).hyp(jj).assignment==tid(cc));
                CL(cl).scan(kk).hyp(jj).assignment(ii)=[]; 
            end
            ii=find(CL(cl).z_id==zid(cc));
            CL(cl).z_id(ii)=[];
            
            CL(ncl+cc).z_id=zid(cc);
            CL(ncl+cc).scan(kk).hyp=strhyp(tid(cc),1);
            CL(ncl+cc).scan(kk).hyp.prob=1;
        end
    end
    end
end





%{

if k-k0 > MHT.nscanback 
        for cl=1:ncl 
            % TBD this should be a while loop for multiple splits
            if ~isempty(CL(cl).scan(k-MHT.nscanback).hyp)
                [split tid2check zidx]=checksplittability(CL,cl,k, MHT);
                if split
                    CL=splitclusters(CL, cl, tid2check, zidx, MHT);
                end
            end
        end
    end
    ncl=size(CL,2);
    
    
    function [split tid2check zidx]=checksplittability(CL, cl, k, MHT)
split=0;
tid2check=0;
zidx=0;
hypmat=cat(1,CL(cl).scan(k-MHT.nscanback).hyp.assignment);

if size(hypmat,2) > 1
    for zz=1:size(hypmat,2)
        tids1=unique(hypmat(:,zz));
        % possibly okay to split with a false measurement not
        % sure!!!
%         tids1=tids1(tids1~=0);
    
        if numel(tids1)==1
            tid2check=tids1;
            split=1;
            for kk=k-MHT.nscanback+1:k-1
                hypmat=cat(1,CL(cl).scan(kk).hyp.assignment);
                idxtid=find(hypmat==tid2check);
                if ~isempty(idxtid)
                [a b]=ind2sub(size(hypmat), idxtid);
                if numel(unique(hypmat(:,b))) > 1
                    split=0;
                end
                if kk==k-1
                    zidx=b(1);
                end
                else
                    split=0;
                end
            end
        end
        % one measurement split only for now....
        if split, break; end
    end
end

function CL = splitclusters(CL, id, tid, zidx, MHT)
% hypotheses
ncl=size(CL,2);
CL(ncl+1)=strClstr(1);
for k=1:size(CL(id).scan,2)
    for jj=1:size(CL(id).scan(k).hyp,2)
        as=CL(id).scan(k).hyp(jj).assignment;
        CL(id).scan(k).hyp(jj).assignment=as(as~=tid);
    end
end

km1=size(CL(id).scan,2);
for k=km1-MHT.nscanback:km1
    CL(ncl+1).scan(k).hyp=strhyp(tid,1);
    CL(ncl+1).scan(k).hyp.prob=1;
end


% current set of measurements
CL(ncl+1).z_id=CL(id).z_id(zidx);
CL(id).z_id=CL(id).z_id(1:end~=zidx);
%}