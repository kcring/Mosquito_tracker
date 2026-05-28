function [tXh, tP, tXh_, tP_, Xh]=shiftdown(k, scan, tXh, tP, tXh_, tP_, Xh, MHT, Xi)
global gf;
global idl;


% ---------------------------------------
% update tracks based on current assignment
% remove all others
% ---------------------------------------
if k > MHT.nscanback
    tk1=cat(1,scan(k-MHT.nscanback).hyp.assignment);
    tk1=tk1(tk1~=0);
    for id=tk1
        [t, br]=id2tbr(id, gf);
        [r, c]=getind(Xi.nX, 1, t, 1:Xi.nX, Xi.nX);
        tXh(r,1,1:end~=br)=0;
        tP(r,c, 1:end~=br)=0;
        idl(t,1:end~=br)=0;
    end
    Xh(:,k-MHT.nscanback)=sum(tXh(:,1,:),3);
    
    % shift
    tXh(:,1:MHT.nscanback,:)=tXh(:,2:MHT.nscanback+1,:);
    tP(:,1:(MHT.nscanback)*Xi.nX,:)=tP(:,Xi.nX+1:(MHT.nscanback+1)*Xi.nX,:);
    
    tXh(:,MHT.nscanback+1,:)=0;
    tP(:,(MHT.nscanback)*Xi.nX+1:(MHT.nscanback+1)*Xi.nX,:)=0;
    
    tXh_(:,1:MHT.nscanback,:)=tXh_(:,2:MHT.nscanback+1,:);
    tP_(:,1:(MHT.nscanback)*Xi.nX,:)=tP_(:,Xi.nX+1:(MHT.nscanback+1)*Xi.nX,:);
    
    tXh_(:,MHT.nscanback+1,:)=0;
    tP_(:,(MHT.nscanback)*Xi.nX+1:(MHT.nscanback+1)*Xi.nX,:)=0;
else
    Xh=[];
end
