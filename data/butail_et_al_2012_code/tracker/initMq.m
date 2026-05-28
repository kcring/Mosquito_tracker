function [tp, Xh]=initMq(Z0,cams, optTrack)
% function [tp T2]=initialize_target(Xh, Z0, Zia, sigma_Z, optTrack)
% 
Xi=strXi;

r=lsTriangulate([Z0(1).u, Z0(2).u], cams);
if r(3) < optTrack.auto.swarm_boundaries(2) && ...
   r(3) > optTrack.auto.swarm_boundaries(1)
    tp=sample_target_state(r, Xi, optTrack);
    [tp wts]=pf_update(Z0, tp, cams, Xi, optTrack);

    tp=mq_motion(tp, optTrack.auto.mm, optTrack.auto.dt);
    Xh=postest(tp,wts, optTrack.auto.pfout);
    
    % perform kmeans clustering to partition the points into the two
    % velocity directions
    %{
    try
    idxclstr=kmeans(tp(Xi.rdi, :)', 2, 'emptyaction', 'singleton');
    catch
        keyboard
    end
    tp1=tp(:,idxclstr==1); 
    tp2=tp(:,idxclstr==2);
    
    % add remaining samples randomly from existing samples
    stp1=size(tp1,2);
    for jj=stp1+1:optTrack.auto.Np
        tp1(:,jj)=tp1(:,(randperm(stp1)==1));
    end
    stp2=size(tp2,2);
    for jj=stp2+1:optTrack.auto.Np
        tp2(:,jj)=tp2(:,(randperm(stp2)==1));
    end
    
    
    % motion model since the target is initialized in the next frame
    tp1=mq_motion(tp1, optTrack);
    tp2=mq_motion(tp2, optTrack);

    Xh=postest(tp1,wts, 1); Xh(Xi.ur)=1;
    Xh2=postest(tp2,wts, 1); Xh2(Xi.ur)=2;
    
    % flag denoting that update was done
    T1=tp1; 
    T2=tp2;
    %}
else
    tp=zeros(Xi.nX, optTrack.auto.Np); 
%     T2=tp;
    Xh=zeros(Xi.nX,1); 
%     Xh2=Xh;
end



function tp = sample_target_state(r, Xi, optTrack)

% Sample particles with high variance in velocity
tp=zeros(Xi.cX(1), optTrack.auto.Np);

switch optTrack.auto.mm
    case 'cv'
        tp(Xi.ri, :)=r*ones(1,optTrack.auto.Np);
        tp(Xi.rdi,:)=randn(3,optTrack.auto.Np)*500;
    case 'wna'
        tp(Xi.ri, :)=r*ones(1,optTrack.auto.Np);
        tp(Xi.rdi,:)=randn(3,optTrack.auto.Np)*500;
        tp(Xi.rddi,:)=randn(3,optTrack.auto.Np)*3000;
    case 'sp1'
        tp(Xi.ri, :)=r*ones(1,optTrack.auto.Np);
        tp(Xi.rdi,:)=normv(randn(3,optTrack.auto.Np));
        tp(7:8,:)=randn(2,optTrack.auto.Np);
        tp(9,:)=randn(1,optTrack.auto.Np)*3000;
    case 'sp2'
        tp(Xi.ri, :)=r*ones(1,optTrack.auto.Np);
        for ii=1:optTrack.auto.Np
            se3rand=orth(rand(3));
            tp(Xi.rdi,ii)=se3rand(:,1);
            tp(Xi.yi,ii)=se3rand(:,2);
            tp(Xi.zi,ii)=se3rand(:,3);
        end
        tp(Xi.si,:)=randn(1,optTrack.auto.Np)*1000;
        tp(Xi.ui,:)=randn(2,optTrack.auto.Np);        
    case 'wienera'
        tp(Xi.ri, :)=r*ones(1,optTrack.auto.Np);
        tp(Xi.rdi,:)=randn(3,optTrack.auto.Np)*500;
        tp(Xi.rddi,:)=randn(3,optTrack.auto.Np)*7000;
    otherwise
        fprintf('[!] Check motion model\n');
end