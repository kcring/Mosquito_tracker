function [Z3d, assoc_t]=mqOccResolve(occ, p, Z3d, assoc_t, cams)

global gdebug;

nc=size(cams,2);
Xi=strXi;

clr=rand(20,3);
markers=['x'; 'o'; 's'; 'd'; 'd';'+'; '*'; '^'; 'x'; 'o'; 's'; 'd'; 'd';'+'; '*';];
mXh=mean(p,2);
dist_t=3;

for oo =1:size(occ,2)
    t=occ(oo).t;
    
    nt=numel(t);
    zh=zeros(2,nt);
    vh=zh;
    for cc=1:nc
        
        silp=unique(cat(1,Z3d(assoc_t(t),cc).pixel_list), 'rows')';
        for t_ii=1:nt
            rt=snipx(Xi.cX, t(t_ii));
            zh(:,t_ii)=w2cam(mXh(rt(Xi.ri),1), cams(cc));
            vh(:,t_ii)=abs(w2cam_rdot(mXh(rt,1), cams(cc), 1/40, Xi));
        end
        
        % call optimization to split the occluded blob into streaks
        x=silp(2,:);
        y=silp(1,:);
        
        r=[x',y'];
        pts0=r';
        model.mu=[];
        for jj=1:nt
            model.mu(jj,:)=zh(:,jj)';
            model.Sigma(:,:,jj)=diag(vh(:,jj));
        end
        
        model.mu=model.mu';
        model.weight=ones(1,nt)/nt;    
        [idx model llh]=emgm(r', model);
        
        
        for ff=1:nt
            pts=[r(idx==ff,1)'; r(idx==ff,2)'];
            
            % Soft cluster
            npts=size(pts,2);
            if npts>3
                A=[pts(1,:)', ones(npts,1)]; b=pts(2,:)';
                l=(A'*A)\A'*b;
                dist=abs(-l(1)*pts0(1,:)+pts0(2,:)-l(2))/sqrt(l(1)^2+1);

                idx1=find(dist<dist_t);


                if gdebug
                    figure(cc+2); gcf;
                    plot(pts0(1,idx1), pts0(2,idx1), markers(ff), 'color', clr(ff,:), ...
                                'markersize',4);
                    hold on;
                    xr=[min(pts0(1,:)) max(pts0(1,:))];                
                    plot(xr, l(1)*xr+l(2), 'g');
                end
            end
        end
    end
 
end

if gdebug
    keyboard
end


