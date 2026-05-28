function hyp=get_newid(nt)
global gf;
global idl;

avids=find(idl(:,1)==0);
aid=avids(1:nt);
tid=tbr2id(aid,1,gf);
hyp=strhyp(tid',1);
idl(aid,1)=1;

hyp.prob=1;