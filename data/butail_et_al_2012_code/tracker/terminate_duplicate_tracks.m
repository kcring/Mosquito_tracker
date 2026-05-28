function [Xh, p, ntem]=terminate_duplicate_tracks(Xh, p, k, optTrack)

% for each option above we look back from the current time-step and see
% if the tracks are "close enough". If they are, we remove those entries
% (in the case of temporary track, the complete track is removed), in case
% of parmanent, only the part that was common is removed

Xi=strXi;
ntem=0;
thresh=100;


ids=getids(Xh,k,Xi);

if numel(ids)>1
    dim=max(ids);

    dist=ones(dim)*100;
    
    for ii=ids'
        for jj=ids(find(ids==ii)+1:end)'
            dist(ii,jj)=trackDist(ii, jj, Xh, k, Xi, optTrack.auto.dup_t);
        end
    end

    idx=find(dist<thresh);
    if ~isempty(idx)
        [t_del1 t_del2]=ind2sub([dim, dim], idx);
        for ii=t_del1'
            % find the shorter track of the two
            r(1,1:Xi.nX)=getind(Xi.nX, k, ii, 1:Xi.nX,1);
            r(2,1:Xi.nX)=getind(Xi.nX, k, t_del2((t_del1==ii)), 1:Xi.nX,1);
            
            tl(1)=sum(Xh(r(1,1),:)~=0);
            tl(2)=sum(Xh(r(2,1),:)~=0);

            [val idx]=min(tl);
            % terminate full track if it is a temporary one
            if val > optTrack.auto.tt_tl
                Xh(r(idx,Xi.nX), k:-1:k-optTrack.auto.dup_t)=0;
            else
                Xh(r(idx,Xi.nX),:)=0;
            end
            p(r(idx,Xi.nX),:)=0;
            ntem=ntem+1;
        end
    end
end



function dist=trackDist(id1, id2, Xh, k, Xi, dup_t)
try
diff=Xh(snipx(Xi.cX, id1, Xi.ri),(k-dup_t+1):k)-Xh(snipx(Xi.cX, id2, Xi.ri),(k-dup_t+1):k);
catch
    keyboard
end
% this ensures that targets which are initiated close by are not terminated
diff(1:3,diff(1,:)==0)=100; 
dist=norm(diff);