function hypmat = assign_new_ids(hypmat)
global gf;
global idl;


newidx=hypmat<gf & hypmat > 0;
nids=unique(hypmat(newidx));
num_new_t=numel(nids);
% new ids are those that have no branches 
avids=find(sum(idl,2)==0);
aid=avids(1:num_new_t);
for ii=1:num_new_t
    hypmat(hypmat==nids(ii))=tbr2id(aid(ii),1,gf);
end
idl(aid,1)=1;
