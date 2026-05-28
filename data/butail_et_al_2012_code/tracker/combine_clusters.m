function CL = combine_clusters(k, CL, Xh_, P_, Zk, MHT, Kalman, Xi)
global gf;
global idl;
nz=size(Zk,1);
ncl=size(CL,2);
% the measurements are current only so all previous are removed
for cl=1:ncl, CL(cl).z_id=[]; end
clgate=ones(nz,ncl)*MHT.gate+1;


for cl =1:ncl
     tids=unique(cat(1,CL(cl).scan(k-1).hyp.assignment)); % all hypothesized targets
     tids=tids(tids~=0);
     CL(cl).kgate=zeros(nz, numel(tids));
     for z=1:nz


        for t=1:numel(tids)
            [tr br]=id2tbr(tids(t),gf);
            [r c]=getind(Xi.nX, k, tr, 1:Xi.nX, Xi.nX);
            try
            S=Kalman.H*P_(r,c,br)*Kalman.H'+Kalman.R;
            catch
                keyboard
            end

            CL(cl).kgate(z,t)=MHT.gatefun((Xh_(r(Xi.ri),k,br)), Zk(z,:)', S);
        end
        clgate(z,cl)=min(CL(cl).kgate(z,:));

        if clgate(z,cl) < MHT.gate
            CL(cl).z_id=[CL(cl).z_id z];
        end
    end
end

% right now there is no good way to ensure that all sets are covered,
% but let's run it for a significant number of times
for runs=1:3 % this can be as many as possible ...
    for cl1=1:ncl
        for cl2=cl1+1:ncl
            if sum(ismember(CL(cl1).z_id, CL(cl2).z_id))
%                     CLT=CL;
                CL=combinehyp(CL, cl1, cl2, MHT);
%                     try
%                         cat(1,CL(cl1).scan(k-1).hyp.assignment)
%                         catch
%                             keyboard;
%                         end
            end
        end
    end
end
ezidx=[];
for cl=1:ncl, if isempty(CL(cl).z_id), ezidx=[ezidx, cl]; end, end
CL(ezidx)=[];


function CL = combinehyp(CL, id1, id2, MHT)

CL(id1).z_id=unique([CL(id1).z_id, CL(id2).z_id]);
km1=size(CL(id1).scan,2);
for k=km1:-1:1
    nh1=size(CL(id1).scan(k).hyp,2);
    nh2=size(CL(id2).scan(k).hyp,2);
    hypnew=strhyp(ones(nh1*nh2,1), 1);
    jj=1;
    if ~isempty(CL(id2).scan(k))
        for h1=1:nh1;
            for h2=1:nh2
                hyp1=CL(id1).scan(k).hyp(h1);
                hyp2=CL(id2).scan(k).hyp(h2);
                hypnew(jj)=combinehypotheses(hyp1, hyp2);
                jj=jj+1;
            end
        end
    end
    if nh1*nh2>=1
        CL(id1).scan(k).hyp=hypnew;
    end
end


CL(id2).z_id=[];
CL(id2).scan=[];

function hyp=combinehypotheses(hyp1, hyp2)

hyp=strhyp([hyp1.assignment hyp2.assignment], ...
                hyp1.parent); %[hyp1.parent, hyp2.parent]
hyp.prob=hyp1.prob*hyp2.prob;
