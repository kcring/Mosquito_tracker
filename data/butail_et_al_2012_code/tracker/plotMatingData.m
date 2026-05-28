function plotMatingData

if nargin < 1
    [floc, frmloc]=getexpdirs;
end

[optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams]=initproc(floc, frmloc);

swarm=1;
mode=3;
datfile1 = [dataloc, sprintf('data_mq_%s.mat', expname)];
if mode==3
    datfile1 = [dataloc, sprintf('data_mq_%s_F_W.mat', expname)];
end

load(datfile1);
metalist=dir([floc, 'calib/metadata_*.txt']);
if(size(metalist,1))

    metadata=csvread([floc, 'calib/', metalist(1).name]);
    if numel(metadata) ==5
        focal_male=metadata(5);
        female=metadata(4);
    else
        error('male female ids not found !');
    end
else
    error('metadata not found !');
end

% you can change camid to 1 or 2
camid=1;
record=1;
if record, visible=0; else visible=1; end

% ---- non-zero time-steps
nz=find(sum(Xh)~=0);
k0=nz(1); kF=nz(end);

     
warning('off', 'Images:initSize:adjustingMag');
cmap=rand(size(Xh,1)/Xi.nX, 3);

r1=Xh(1:Xi.nX:end,:);
r2=Xh(2:Xi.nX:end,:);
r3=Xh(3:Xi.nX:end,:);


if mode==3, factor=1000; else factor=1; end

Xh(Xh==0)=nan;
alim=[min(r1(r1~=0)) max(r1(r1~=0)) min(r2(r2~=0)) max(r2(r2~=0)) min(r3(r3~=0)) max(r3(r3~=0))]/1000;
for k=k0:kF
    if visible
        figure(1); gcf; clf;
    else
        figure('Visible', 'off'); clf;
    end
    if mode==2, img=imread([frmloc, frmlist(k,camid).name]); end
    tail=max(k0, k-5);
    nz_mq=find(~isnan(Xh(1:Xi.nX:end,k)));
    
    if ~swarm, nz_mq=intersect(nz_mq, [focal_male, female]'); end
    if mode==2 && ~isempty(nz_mq)
        try
        [img trks colors]=highlightMq(img, nz_mq, tail, k, Xh, cams(camid), Xi, cmap);
        catch
            keyboard
        end
        imshow(img); hold on; 
    end
    
    for jj=nz_mq'
        
        r=getind(Xi.nX, k, jj, Xi.ri, 1);
        if mode==3
            colr=cmap(jj,:)+.5*(ones(1,3)-cmap(jj,:));
            if jj==focal_male, colr=[0 0 1]; end
            if jj==female, colr=[1 0 0]; end
        else
            colr=colors(jj,:);
            if jj==focal_male, colr=[0 0 1]; end
            if jj==female, colr=[1 0 0]; end
        end
        if mode==3
            trk=Xh(r,tail:k)/factor;
            plot3(trk(1,:), trk(2,:), trk(3,:), 'Color', colr, 'linewidth', 2);
            hold on;
            plot3(trk(1,end), trk(2,end), trk(3,end), 'o', 'Color', colr, 'linewidth', 2);
        else
            plot(trks(2*jj-1,:), trks(2*jj,:), 'Color', colr);
            plot(trks(2*jj-1,end), trks(2*jj,end), 'Color', colr);
            
%             text(trks(2*jj-1,end)+3, trks(2*jj,end)+3, sprintf('%d', jj), 'Color', colors(jj,:));
        end
    end
    if mode==3, 
        axis(alim); 
        set(gca, 'fontsize', 16);
        box on;
    end
    if ~record, title(k); end
    drawnow;
    if(record)
        set(gcf,'units','pixel');
        set(gcf,'position',[0,0,1392/1.667,1040/1.667]);
        set(gcf,'papersize',[1392/1.667,1040/1.667]);
        set(gcf,'PaperPositionMode','auto')
        %                 set(gca,'LooseInset',get(gca,'TightInset'));
        if mode==3
            set(gca, 'OuterPosition', [0 0 1 1]);
        else
            set(gca, 'OuterPosition', [-.18 -.18 1.35 1.35]);
        end
        print('-dbmp', sprintf('/tmp/mating_cam%d_%.4d.bmp', camid, k));
    end
    if ~record, pause(.01); end
end
    

function [img trks colors]=highlightMq(img, nz_mq, tail, k, Xh, cam, Xi, cmap)

[h, w]=size(img(:,:,1));

% set the circular region for each mosquito
npts=200; th=linspace(0,2*pi, npts);
hl_rad=0:15;
[HLR TH]=meshgrid(hl_rad, th);
regionx0=ceil(HLR.*cos(TH));
regiony0=ceil(HLR.*sin(TH));


colors=zeros(max(nz_mq),3);
for jj=nz_mq'
    r=getind(Xi.nX, k, jj, Xi.ri, 1);
    trks(2*jj-1:2*jj,:)=w2cam(Xh(r,tail:k), cam);
    colors(jj,:)=cmap(jj,:)+.5*(ones(1,3)-cmap(jj,:));
end
for jj=1:size(trks,1)/2
    region_ctr=[trks(2*(jj-1)+1,end), trks(2*jj,end)]';
    regionx=ceil(region_ctr(1))+regionx0;
    regiony=ceil(region_ctr(2))+regiony0;
    hlght=unique([regionx(:), regiony(:)], 'rows');

    hlght(:,2)=min(h,hlght(:,2));
    hlght(:,2)=max(1,hlght(:,2));
    hlght(:,1)=min(w,hlght(:,1));
    hlght(:,1)=max(1,hlght(:,1));
    % brighten up this part in the image
    imgind=sub2ind([h,w], hlght(:,2), hlght(:,1));
    [num thresh]=hist(double(img(imgind)),5);
    thresh=mean(thresh);
    img(imgind)=img(imgind)*1.1;
    for gg=1:numel(imgind)
        if img(imgind(gg))<thresh(1)
            img(imgind(gg))=img(imgind(gg))*.99;
        end
    end
end

