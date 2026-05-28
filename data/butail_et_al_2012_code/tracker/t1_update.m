function [Xh instr]= t1_update(Xh, mqid, kt, frame, get_cam_calib, Xi, Zi, optTrack)

Pr=[];
jj=mqid;

% Cylinder points 
% To be changed to a truncated cone later
rc=0:50;
tc=linspace(-pi, pi, 30);

startp=optTrack.trackone.swarm_d0;
endp=optTrack.trackone.swarm_d0+optTrack.trackone.swarm_s0;


[r c]=getind(Xi.nX, kt, mqid, Xi.ri, 1);
if(size(Xh,1)>=r(1))
    if(Xh(r(1),c))
    %     startp=norm(frame(k-1).mq(jj).X(1:3))-300;
        startp=norm(Xh(r,c))-300;
        endp=startp+500;
    end
end

switch optTrack.trackone.alg
    case 'mpf'
        dc=0:5:(endp-startp);
        [Rc Tc Dc]=meshgrid(rc,tc,dc);

        Xc=Rc.*cos(Tc);
        Yc=Rc.*sin(Tc);
        Zc=Dc;

        %%% projecting epipolar line %%

        lcam=get_cam_calib(1);
        rcam=get_cam_calib(2);

        % store mosquito values in a local variable
        if(isfield(handles, 'frame'))
            if(size(frame, 2)>=k)
                lc=frame(k).mq(jj).c(:,1);
            end
        end

        axes(ah(1));

        if(exist('lc', 'var'))
            eline_w=get_eline([lc(1), lc(2)], ...
                        lcam, startp, endp);


            eline_rcam=w2cam(eline_w, rcam);   
%             frame(k).reline(:,:,jj)=eline_rcam;
            % create a region where you'll look for mosquitoes in the right camera
            % frame

            % slope is (y2-y1)/(x2-x1)
            m0=(eline_rcam(2,2)-eline_rcam(2,1))/...
               (eline_rcam(1,2)-eline_rcam(1,1));
            mp=-1/m0;

            % intercepts
            % y=mx+c; at x=0 y=c, therefore c-y1/(0-x1)=m;
            c0=-eline_rcam(1,1)*m0+eline_rcam(2,1);
            cp=-eline_rcam(1,1)*mp+eline_rcam(2,1);

            m1=m0;m3=m0;
            m2=mp;m4=mp;
            c1=c0-optTrack.bbox(1);
            c3=c0+optTrack.bbox(1);

            % swarm size in pixels
            ssp=sqrt((eline_rcam(1,end)-eline_rcam(1,1))^2 + ...
                     (eline_rcam(2,end)-eline_rcam(2,1))^2);
            c2=cp;
            c4=-eline_rcam(1,end)*mp + eline_rcam(2,end);

            % number of dots in right cam
            rcam_dots=c(2).centroids;

            nce=size(rcam_dots,1);
            if(optTrack.debug)
                axes(ah(2));
                plot(rcam_dots(:,1), rcam_dots(:,2), 'r.', 'MarkerSize', ms(2));
            end

            eline_v3=(eline_w(:,2)-eline_w(:,1))/norm(eline_w(:,2)-eline_w(:,1));
            % making an orthogonal frame since the directions of other two don't
            % matter we choose v1=[a b c], v2=[0 -c -b] and v3= v1xv2;
            eline_v2=[0 -eline_v3(3) eline_v3(2)]'/norm([0 -eline_v3(3) eline_v3(2)]);
            eline_v1=cross(eline_v2, eline_v3);
            wTeline=[eline_v1, eline_v2, eline_v3, eline_w(:,1);
                        0 0 0 1];

            % points in 3D where we need to find probability values
            [re1 re2 re3]=trpa2b(Xc, Yc, Zc, wTeline);
            [uCL vCL]=wp2cam(re1, re2, re3, lcam);
            Pr_camL=normpdf(uCL, lc(1), optTrack.trackone.camstd(1)).*...
                    normpdf(vCL, lc(2), optTrack.trackone.camstd(1));

            [uCR vCR]=wp2cam(re1, re2, re3, rcam);

            if(optTrack.debug) % set this value to 1 if you wish to see how it all looks
                axes(ah(2));
                xt=400:900;
                plot(xt, m1*xt+c1, '--b');
                plot(xt, m2*xt+c2, '--r');
                plot(xt, m3*xt+c3, '--g');
                plot(xt, m4*xt+c4, '--y');
            end

            zz=0;
            for ii=1:nce
                xi=rcam_dots(ii,1);
                yi=rcam_dots(ii,2);

                % to fit in a tilted box we ask that the y values be between the
                % enclosing lines m1 is the line parallel to epipolar line just
                % below it (by optTrack.bbox) m2 is normal to it at one end, m3 is
                % again parallel and above the epipolar line, and m4 is normal at
                % the other end

                % also since we know that for any 3D point visible in left frame
                % the same point will have a lesser x-value in the right frame, we
                % place that condition as well

                if(m1>0 && yi>m1*xi+c1 && ...
                        yi>m2*xi+c2 && ...
                        yi<m3*xi+c3 && ...
                        yi<m4*xi+c4 || ...
                   m1<0 && yi>m1*xi+c1 && ...
                        yi<m2*xi+c2 && ...
                        yi<m3*xi+c3 && ...
                        yi>m4*xi+c4)
                    % additional condition that the mosquito pixel diff as per the
                    % baseline change more than 30 pixels in each x-y direction.
                    % The reasoning for this is that a mosquito may not move too
                    % fast in the direction of the optical axis. Thus the pixel
                    % difference f*b/z is not more than 30 pixels in the x
                    % direction and 10 pixels in the y direction
                    zz=zz+1;
                    rcam_selected_dots(:,zz)=[xi, yi]';
                end
            end

            for ii=1:zz
                xi=rcam_selected_dots(1,ii);
                yi=rcam_selected_dots(2,ii);
                Pr_camR=normpdf(uCR, xi, optTrack.trackone.camstd(2)).*...
                    normpdf(vCR, yi, optTrack.trackone.camstd(2));
                pdf=Pr_camL.*Pr_camR;
                [valm idxm]=max(pdf(:));
                [i1 i2 i3]= ind2sub(size(pdf), idxm);
                map_r(:,ii)=[re1(i1,i2,i3), re2(i1,i2,i3), re3(i1,i2,i3)]';
                Pr(ii)=sum(sum(sum(pdf)));    
            end

            frame(k).mq(jj).eline=eline_rcam;
            plot_eline(handles, 2);

            if(~isempty(Pr))
                [val idx]=sort(Pr);
                num_Pr=length(Pr);

                % choose top 5 or num_Pr whichever is less
                if (num_Pr <5)
                    frame(k).mq(jj).cam(2).candpos=rcam_selected_dots(:,idx);
                else
                    frame(k).mq(jj).cam(2).candpos=rcam_selected_dots(:,idx(end-4:end));
                end

                    frame(k).mq(jj).c(:,2)=rcam_selected_dots(:,idx(end));

                    [r c]=getind(Zi.nZ, k, mqid, Zi.uv, 2);
                    Z(r,c(2))=rcam_selected_dots(:,idx(end));
                    [r c]=getind(Zi.nZ, k, mqid, Zi.ua, 2);
                    Z(r,c(2))=1;

                    [r c]=getind(Xi.nX, k, jj, 1:Xi.nX, 1);
                    Xh(r(Xi.ri), c)=map_r(:,idx(end));

                    % *** if backtracking 
                    if (backtrack)
                        if(Xh(r(1), c+1))
                            Xh(r(Xi.rdi),c)=(Xh(r(Xi.ri),c+1)-Xh(r(Xi.ri),c))/optTrack.dt;
                         end
                    else
                        if(Xh(r(1), c-1))
                            Xh(r(Xi.rdi),c)=(Xh(r(Xi.ri),c)-Xh(r(Xi.ri),c-1))/optTrack.dt;
                        end
                    end
                    % nf(2) == 2 for manual tracking
                    Xh(r(Xi.fi), c)=[-1; 2];
                    instr=sprintf('3D position: [%d, %d, %d]\nPress next (>) or Record to store all position information', ...
                                    ceil(map_r(1,idx(end))), ceil(map_r(2,idx(end))), ceil(map_r(3,idx(end))));
                    

            else
                instr='[!] No match found in the right camera image... Mark a point manually (check binary threshold)';
                % 
            end
        end
end



function pts = get_eline(uv, cam1, startp, endp)

% cam1 is the frame where the points are located
% cam2 is the frame where lines are to be drawn

% get distorted values in mm in camera frame
dist=inv(cam1.km)*[uv(1) uv(2) 1]';

% convert to world coordinates

% set optimization options
options = optimset('Display','off','TolFun',1e-6);
[xyn, fval, residual]=fmincon(@(xyn) undist(xyn, dist(1:2), cam1.kc1, cam1.kc2), ...
                        dist(1:2), ... 
                        [],[],[],[],...
                        [-2, -2],[2, 2,],...
                        [], options);
% points from 0 along this vector to the end of the tank

% normalize
xyn=[xyn;1];
xyn=xyn/norm(xyn);
pts=xyn*(startp:1:endp);

% put it in world coordinates
pts=tra2b(pts, inv(cam1.trm));
