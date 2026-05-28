function p=update_and_initialize(k, zt, t, Z3d, p, p_, PF, Xi)
global gf;
[t br]=id2tbr(t, gf);
[r c]=getind(Xi.nX, k, t, 1:Xi.nX, PF.Np);
if size(p_,1)>=r(1) % if this is an existing target
    tree=br;
    for z=zt'
        if p_(r(1),c(1),tree)
            p(r,c, br)= PF.lfn(Z3d(z,:), p_(r,c, tree));
        end
        br=br+1;    
    end
else
    % if its a new target it should not be assigned to multiple
    % measurements
    p(r,c, 1)=PF.inittarget(Z3d(zt,:));
end
