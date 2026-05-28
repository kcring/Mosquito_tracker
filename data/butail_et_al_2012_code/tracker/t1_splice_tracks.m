function [Xh, splice_flag] = t1_splice_tracks(Xh, Xh0, id, k, dir, Xi, Xi0, thresh)
% this function takes in a data matrix and attaches the track that is
% closest to it within a threshold
splice_flag=0;
% get the 3D position in the current time step
[r, c]=getind(Xi.nX, k, id, 1:Xi.nX, 1);
Xh_id=Xh(r,c);
% kr1=find(Xh(snipx(Xi.cX, id, Xi.ri(3)),:)~=0);
% get the 3D positions of all automatically generated tracks in the current
% time step

mqauto_k=show_tracked_mq(Xh0, Xi0, k, 0);


% compare distance (and velocity) within threshold
jj=1;
dist=zeros(numel(mqauto_k),1);
for mq_ii=mqauto_k'
    [r0, c0]=getind(Xi0.nX, k, mq_ii, 1:Xi0.nX, 1);
    dist(jj)=norm(Xh_id(1:3)-Xh0(r0(Xi.ri), c0));
    jj=jj+1;
end

[val, idx]=min(dist);
% fprintf('[I] Closest track is %.2f mm away...\n', val);

if val < thresh
    if ~dir
        r0=getind(Xi0.nX, k, mqauto_k(idx), 1:Xi0.nX, 1);
        kr2=find(Xh0(r0(1),:)~=0);
        % change splicing should take place from the next frame onwards
        Xh(r(Xi.ri),k+1:kr2(end))=Xh0(r0(Xi.ri),k+1:kr2(end));
        Xh(r(Xi.rdi),k+1:kr2(end))=Xh0(r0(Xi.rdi),k+1:kr2(end));
%         Xh(snipx(Xi.cX,id,Xi.fi(2)),kr1(1):kr2(end))=2;
%     else % backwards suggestion for splicing not to be done
%         kr2=find(Xh0(snipx(Xi.cX, mqauto_k(idx), Xi.ri(3)),:)~=0);
%         Xh(snipx(Xi.cX,id),kr2(1):k)=Xh0(snipx(Xi.cX, mqauto_k(idx)),kr2(1):k);
%         Xh(snipx(Xi.cX,id,Xi.fi(2)),kr2(1):kr1(end))=2;
    end
    splice_flag=1;        
end

