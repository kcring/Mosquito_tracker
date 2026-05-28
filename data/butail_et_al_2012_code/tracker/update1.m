function [tXh tP hypmat]=update1(k, hypmat, zt, t, Zk, tXh, tP, tXh_, tP_, MHT, Kalman, Xi)
global gf;
global idl;
k1=min(k, MHT.nscanback+1); % to search the temporary state space

t0=t;
[t br]=id2tbr(t, gf);
[r c]=getind(Xi.nX, k1, t, 1:Xi.nX, Xi.nX);

tree=br;
for z=zt'
    if tXh_(r(1),k1,tree)
        [tXh(r,k1,br), tP(r,c, br)]= ...
                    kalmanUpdate(tXh_(r,k1,tree), tP_(r,c, tree), Zk(z,:)', Kalman);
    end
    try
    tid=tbr2id(t,br,gf);
    idl(t,br)=1;
    hypmat(hypmat(:,z)==t0,z)=tid;
    abr=find(idl(t,:)==0);
    br=abr(1);  
    
    catch
        keyboard
    end
end

