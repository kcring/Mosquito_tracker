function Xh= cleanup_space(k, Xh, Xi)
global gf;
global idl;

tl=sum(Xh(1:Xi.nX:end,:)~=0,2);

if ~isempty(tl)
    tle=Xh(1:Xi.nX:end,end)~=0;
    idno=find(tl==1 & tle==0);
    for ii=idno'
        idl(ii,:)=0;
        r=getind(Xi.nX, 1, ii, 1:Xi.nX,1);
        Xh(r,:)=0;
    end
end