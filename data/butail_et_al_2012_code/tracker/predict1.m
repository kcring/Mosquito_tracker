function p_=predict1(k, p, p_, PF, Xi)

for ll=1:size(p,3)
    for t=1:size(p,1)/Xi.nX
        [r, c]=getind(Xi.nX, k, t, 1:Xi.nX, PF.Np);
        if p(r(1),c(1)-1,ll)
            p_(r,c,ll)=PF.motion(p(r, c-PF.Np, ll));
        end
    end
end

