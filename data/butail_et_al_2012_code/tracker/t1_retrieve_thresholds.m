function [binary_t area_t]=t1_retrieve_thresholds(Z, k, mqid, camid, Zi)

[r c]=getind(Zi.nZ, k, mqid, Zi.bt, 2);
if size(Z,1)>=r && size(Z,2)>=c(2)
    if Z(r,c(camid))
        binary_t=Z(r,c(camid));
    else
        binary_t=0;
    end

    [r c]=getind(Zi.nZ, k, mqid, Zi.at, 2);
    if Z(r,c(camid))
        area_t=Z(r,c(camid));
    else
        area_t=0;
    end
else
    binary_t=0;
    area_t=0;
end

