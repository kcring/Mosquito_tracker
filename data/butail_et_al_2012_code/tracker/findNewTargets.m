function [Xh, p, new_t]= findNewTargets(Xh, p, k, cams, Z3d, unassigned, optTrack)

Xi=strXi;
new_t=0;
% Xh_tii=zeros(Xi.nX,2);
% p_tii=zeros(Xi.nX,optTrack.auto.Np,2);

for ii=unassigned
   Z0=[Z3d(ii,1), Z3d(ii,2)];
   
   % Xh_tii is a double sized vector
%    [p_tii(:,:,1), p_tii(:,:,2), Xh_tii(:,1), Xh_tii(:,2)]=initMq(Z0, cams, optTrack);
    [p_tii, Xh_tii]=initMq(Z0, cams, optTrack);
   
   
   if p_tii(1,1,1)
       jj=1;
%        for jj=1:2
            avid=available_ids(Xh, Xi, optTrack.auto.max_tt);
            % note that the timestep is the maximum because the track in a
            % temporary target goes from current timestep to previous
            [r c]=getind(Xi.cX(1), k, avid, 1:Xi.cX(1), 1);

            % the flags are not set yet. Set them when moving to confirmed
            % targets
            Xh(r,c)=Xh_tii(:,jj);
            p(r,:)=p_tii(:,:,jj);
            % NOTE:if you want two targets per new target then call the avid
            % again and set the p(avid).p=T2
            new_t=new_t+1;
%        end
   end
end