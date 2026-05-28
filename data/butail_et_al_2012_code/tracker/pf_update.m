function [p_, wts]=pf_update(Z, p_, cams, Xi,  optTrack)

wts=ones(1, optTrack.auto.Np);
for cc=1:size(cams,2)
    wts=wts.*p_lfn(Z(cc), p_, cams(cc), Xi, optTrack);
end
% RBPF style approx
Uk=cat(2, Z.u);
if Uk(1,1) && Uk(1,2)
    r=lsTriangulate(Uk, cams);
    p_(Xi.ri,:)=r*ones(1,optTrack.auto.Np);
    
    wts=wts/sum(wts);
    neff=1/sum(wts.^2);
    if (neff <=  optTrack.auto.Np/2)
        p_=p_(:,resample(wts));
    end
end