function scan=clustr(k, scan, p_, Z3d, cams, PF, Xi, t_gate)
global gf;
global idl;

% the clusters are the same
scan = copycluster(k-1, k, scan);


ncl=size(scan(k).C,2);
nz=size(Z3d,1);

vm=zeros(nz, ncl);

% just add measurements to each cluster
for zz=1:nz
    for cl=1:ncl
        gatecheck=vclz(k, p_, scan(k).C(cl), Z3d(zz,:), cams, Xi, PF.Np);
        if  gatecheck && gatecheck< t_gate*2
            scan(k).C(cl).z_id=[scan(k).C(cl).z_id zz];
            vm(zz,cl)=1;
        end
    end
end
for zz=1:nz
    if ~sum(vm(zz,:)) % unassociated measurements
        scan(k).C=[scan(k).C strClstr(1)];
        scan(k).C(end).z_id=zz;
        scan(k).C(end).vm=[];
        tid=tbr2id(idl,1,gf);
        scan(k).C(end).hyp=strhyp(tid,1);
        scan(k).C(end).t_id=tid;
        scan(k).C(end).hyp.prob=1;
        idl=idl+1;
    end
end


       
function scan = copycluster(k1, k2, scan)

% only copy t_ids
scan(k2).C=strClstr(size(scan(k1).C,2));

% don't copy z_ids;
for ii =1:size(scan(k2).C,2)
    scan(k2).C(ii).t_id=scan(k1).C(ii).t_id;
end