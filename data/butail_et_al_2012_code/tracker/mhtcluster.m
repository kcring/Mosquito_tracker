function CL=mhtcluster(k,k0, Xh_, P_, CL,Zk, MHT, Kalman, Xi)
if k>k0
    
    CL=combine_clusters(k, CL, Xh_, P_, Zk, MHT, Kalman, Xi);
    
    CL=split_clusters(k, CL, Xh_, P_, Zk, MHT, Kalman, Xi);
    
    CL=add_clusters(CL, Zk);
    
else
    % a cluster per measurement in the first step
    nz=size(Zk,1);
    CL=strClstr(nz);
    for cl=1:nz
        CL(cl).z_id=cl;
    end
end



