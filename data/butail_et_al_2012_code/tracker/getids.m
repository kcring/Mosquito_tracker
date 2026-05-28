function ids=getids(Xh, k, Xi)
ids=find(Xh(Xi.ri(1):Xi.cX(1):end,k)~=0);
if isempty(ids), ids=[]; end