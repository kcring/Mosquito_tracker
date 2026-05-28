function [allwts, Xh, p, unassigned]=mqda_custom(Xh, p, Z3d, cams, k, optTrack)
% function [assoc_t, occ, imstream, wtvals, Xh, p]=mqda_custom(Xh, p, imstream, cams, k, tcolors, optTrack)%
% custom because we use adaptive thresholding to find new measurements and
% reassign them

Xi=strXi;
ids=getids(Xh,k ,Xi);

nc=size(Z3d,2);
allwts=ones(max(ids),size(p,2));

% ------------------------
% Data association
% ------------------------
nZk=size(Z3d,1);
costmat=ones(max(ids), nZk)*9999;
assoc_t=zeros(max(ids),1); 
cost_t=ones(max(ids),1)*9999;

for t_ii=ids'
    r=getind(Xi.nX, k, t_ii, 1:Xi.nX, 1);
%     S=cov(p(r(Xi.ri),:)');
%     mp=mean(p(r(Xi.ri),:),2);
    for zz=1:nZk
        costmat(t_ii,zz)=geticov(p(r,:), Xi, Z3d(zz,1)', cams(1), optTrack) + ...
                        geticov(p(r,:), Xi, Z3d(zz,2)', cams(2), optTrack);
    end
%     geticov(p(r,:), Xi, Z3d(:,1)', cams(1), optTrack) + ...
%                         geticov(p(r,:), Xi, Z3d(:,2)', cams(2), optTrack);
end

ctmp=costmat(ids,:);
switch optTrack.auto.da
    case 'nnda' % nearest neighbor data association
        [cost_t(ids), assoc_t(ids)]=min(ctmp, [], 2);
    case 'gnn' % hungarian method / global nearest neighbor
        a=munkres(ctmp);
        % if a row has all zeros, assoc_t is still 
        [val, assoc_t(ids)]=max(a,[],2);
        assoc_t(ids(sum(a,2)==0))=0;

        % since gnn is one-to-one, the zero assignments automatically
        % cancel out due to high cost initialization
        for ii=ids'
            if assoc_t(ii)
                cost_t(ii)=ctmp(ids==ii,assoc_t(ii));
            end
        end
        
    case 'lgnn' % local gnn
        % note the transpose because clustrgen assumes targets are columns
        ctmp=ctmp';
        
        vm=[ones(nZk,1), genvmat(ctmp, 2*optTrack.auto.gatesize)];
        C=clustrgen(vm); 
        assoc_tmp=zeros(numel(ids),1);
        for cl=1:size(C,2)
            if numel(C(cl).z_id) > 5
                fprintf('nz=%d  ...\n', numel(C(cl).z_id));
            end
            assoc_tmp=lgnnda(C(cl).t_id, C(cl).z_id, ctmp, assoc_tmp, optTrack.auto.gatesize );
        end
       
        % finally get the ids back
        assoc_t(ids)=assoc_tmp;
        
        % since gnn is one-to-one, the zero assignments automatically
        % cancel out due to high cost initialization
        for ii=ids'
            if assoc_t(ii)
                cost_t(ii)=ctmp(assoc_t(ii), ids==ii);
            end
        end
        
    case 'mht0'
        % note the transpose because clustrgen assumes targets are columns
        ctmp=ctmp';
        
        vm=[ones(nZk,1), genvmat(ctmp, 2*optTrack.auto.gatesize)];
        C=clustrgen(vm); 
        assoc_tmp=zeros(numel(ids),1);
        for cl=1:size(C,2)
            if numel(C(cl).z_id) > 7
                assoc_tmp=lgnnda(C(cl).t_id, C(cl).z_id, ctmp, assoc_tmp, optTrack.auto.gatesize );
                fprintf('nz=%d using lgnn ...\n', numel(C(cl).z_id));
            else
                at=assoc_tmp;
                assoc_tmp=mht0da(C(cl).t_id, C(cl).z_id, vm, at, Z3d, cams, p, Xi);
                assoc_tmp1=lgnnda(C(cl).t_id, C(cl).z_id, ctmp, at, optTrack.auto.gatesize );
                if norm(assoc_tmp-assoc_tmp1)
%                     assoc_tmp';
%                     assoc_tmp1';
%                     keyboard
                end
            end
        end
        
        % finally get the ids back
        assoc_t(ids)=assoc_tmp;
        
        % since gnn is one-to-one, the zero assignments automatically
        % cancel out due to high cost initialization
        for ii=ids'
            if assoc_t(ii)
                cost_t(ii)=ctmp(assoc_t(ii), ids==ii);
            end
        end
    case 'mht1'
        
    case 'jpda'
        % probabilistic for each measurement
end


% ------------------------
% Condition to terminate targets for nn, gnn only
% if a target is not associated then kill it
% (with a laser!)
% ------------------------
switch optTrack.auto.da
    case {'nnda', 'gnn'}
        for t_ii=ids'
            % first terminate condition (based on association)
            % may happen if gnn is used
            if cost_t(t_ii) > 2*optTrack.auto.gatesize
                  assoc_t(t_ii)=0;
            end
        end
end


% ------------------------
% Occlusion handling
% ------------------------
% get ids from measurement pairs
if nZk > 0
zpercam=[cat(1, Z3d(:,1).id)'; cat(1, Z3d(:,2).id)'];
aspercam(:,find(assoc_t~=0))=zpercam(:,assoc_t(assoc_t~=0));
% remove temporary ids
for t_ii=ids'
    r=getind(Xi.nX, k, t_ii, 1:Xi.nX, 1);
    tl=sum(Xh(r(1),:)~=0);
    if tl <= optTrack.auto.tt_tl, aspercam(:,t_ii)=0; end
end
% find occlusions
[occ noc]=occlusion_reasoning(aspercam);
if noc
    [Z3d, assoc_t]=mqOccResolve(occ, p, Z3d, assoc_t, cams);
end
end


% nt_terminate=ids(sum(assoc_t(ids),2)==0 & Xh(ids*Xi.ur,k)==0);
% the above method doesn't really work since the mosquitoes move a lot in a
% frame
nt_terminate=ids(sum(assoc_t(ids),2)==0);
mfpr(k,sprintf('%d targets terminated [unassociated]... ', numel(nt_terminate)));
fprintf('terminated target ids=');
fprintf('%d...', nt_terminate);
fprintf('\n');
for t_ii=nt_terminate'
    r=snipx(Xi.cX, t_ii);
    tl=sum(Xh(r(1),:)~=0);
    if tl > optTrack.auto.tt_tl
        Xh(r,k)=0;
    else
        Xh(r,:)=0;
    end
    p(r,:)=0;
end
ids=getids(Xh,k,Xi);


% ------------------------
% Update step
% ------------------------
r1=Xh(1:Xi.nX:end,k); r1=r1(r1~=0);
r2=Xh(2:Xi.nX:end,k); r2=r2(r2~=0);
r3=Xh(3:Xi.nX:end,k); r3=r3(r3~=0);
mpos=mean([r1,r2,r3])';
cpos=cov([r1,r2,r3]);
for t_ii=ids'
    
    r=getind(Xi.nX, 1, t_ii, 1:Xi.nX,1);
    
    if assoc_t(t_ii)
        Z=Z3d(assoc_t(t_ii),:); 
        [p(r,:) wts]=pf_update(Z, p(r,:), cams, Xi, optTrack);

        Xh(r,k)=postest(p(r,:), wts, optTrack.auto.pfout);     
        allwts(t_ii,:)=wts;
    end
end

% -------- Unassigned measurements
unassigned=find(~ismember(1:nZk, assoc_t));


function assoc_tmp=lgnnda(t_id, z_id, ctmp, assoc_tmp, gatesize)

t=t_id; t=t(t~=0);
z=z_id; z=z(z~=0);
if numel(z) > 1
    [a, cost_cl]=munkres(ctmp(z,t));
else
    [cost_cl, a1]=min(ctmp(z,t));
    a=zeros(1,numel(t)); a(a1)=1;
end

if cost_cl > 2*gatesize*numel(t)
    fprintf('[!] Cost of a cluster assignment should not be more than threshold\n');
end
% for each target you'll get a zindex
[val, zidx]=max(a,[],1);
assoc_tmp(t)=z(zidx);
% if a target was not assigned a measurement don't associate
assoc_tmp(t(sum(a,1)==0))=0;


function assoc_tmp=mht0da(t_id, z_id, vm, assoc_tmp, Z3d, cams, p, Xi)

t=t_id; t=t(t~=0);
z=z_id; z=z(z~=0);
vm1=vm(z,t+1);
hyp=hypgen2([ones(size(vm1,1),1), vm1], t);


hfun=@(x,c) (w2cam(x(1:3,:), cams(c)));
prob=zeros(size(hyp,1),1);
for hh=1:size(hyp,1)
    prob(hh)=hypprob(hyp(hh,:), t, Z3d(z,:), p, hfun, Xi);
end
prob=prob/sum(prob);
[val idx]=max(prob);

% measurement-based
for jj=1:size(hyp,2)
    if hyp(idx,jj) % if its a false alarm that target automatically goes
        assoc_tmp(hyp(idx,jj))=z(jj);
    end
end
