function chi2=vclz(k, p_, C, z3d, cams, Xi, Np)
global gf;
jj=1;
chi2=0;

for t=C.t_id
    [tr br]=id2tbr(t,gf);
    [r c]=getind(Xi.nX, k, tr, 1:Xi.nX, Np);
    for v=1:2
        pzh_=w2cam(p_(r(Xi.ri),c, br), cams(v));
        S=cov(pzh_');
        zh_=mean(pzh_,2);
        chi2(jj)=chi2(jj)+(zh_-z3d(v).u)'/S*(zh_-z3d(v).u);

        if det(S) < 10^-5
            %keyboard;
            fprintf('[!] Badly scaled innovation matrix... \n');
            % catching badly scaled matrix
        end    
    end
    jj=jj+1;
end
chi2=min(chi2);