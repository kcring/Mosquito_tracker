function CL=add_clusters(CL, Zk)
nz=size(Zk,1);
ncl=size(CL,2);
zidx=cat(2, CL.z_id);
newz=find(~ismember(1:nz, zidx));
CL(ncl+1:ncl+numel(newz))=strClstr(numel(newz));
for jj=1:numel(newz)
    CL(ncl+jj).z_id=newz(jj);
end