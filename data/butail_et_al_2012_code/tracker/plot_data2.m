function plot_data2
%function plotMqData(datafile, frmloc)
%
% if only datafile is passed then 3D is plotted
% if frmloc is also passed then the estimates are reprojected
% on the screen
% 
% mode=1 plot all data on a 3D plot
% mode=2 write stereo images

datafile = '/Users/puneet2895/Dropbox/marker_exp/data/20160325_191717/output/data/data_mq_EXP_2016_03_25_L_19.17.mat';
frmloc = '/Volumes/Puneet_Back/20160325_191717/frames/';
mode=2; camid=2;
% 
% frmlist = dir(frmloc);
% frmlist = frmlist(end-3000:end-1);
% c1 = get_cam_calib_mar252016(1);
% c2 =  get_cam_calib_mar252016(2);

% cams = 
if mode==2
    load(datafile, 'Xh', 'Xi', 'cams', 'frmlist');
else
    load(datafile, 'Xh', 'Xi');
end
% you can change camid to 1 or 2

record=1;
if record, visible=0; else visible=1; end

% ---- non-zero time-steps
nz=find(sum(Xh)~=0);
k0=nz(1); kF=nz(end);
compass=[.25*cos(0:.1:2*pi);
         .25*sin(0:.1:2*pi)];
     
warning('off', 'Images:initSize:adjustingMag');
cmap=rand(size(Xh,1)/Xi.nX, 3);

r1=Xh(1:Xi.nX:end,:);
r2=Xh(2:Xi.nX:end,:);
r3=Xh(3:Xi.nX:end,:);


if mode==3, factor=1000; else factor=1; end

Xh(Xh==0)=nan;
alim=[min(r1(r1~=0)) max(r1(r1~=0)) min(r2(r2~=0)) max(r2(r2~=0)) min(r3(r3~=0)) max(r3(r3~=0))]/factor;

if mode==1
    r1(r1==0)=nan;
    r2(r2==0)=nan;
    r3(r3==0)=nan;
    figure(1); gcf; clf;
    plot3(r1', r2', r3');
end

if mode>1
    for k=k0:kF
        if visible
            figure(1); gcf; clf;
        else
            figure('Visible', 'off'); clf;
        end
        if mode==2
            img1=imread([frmloc, frmlist(k).name]); 
            img2=imread([frmloc, frmlist(k+1500).name]); 
        end
        tail=max(k0, k-5);
        nz_mq=find(~isnan(Xh(1:Xi.nX:end,k)));


        if mode==2
            subplot(1,2,1);
            [img, trks, colors]=highlightMq(img1, nz_mq, tail, k, Xh, c1, Xi, cmap);
            imshow(img); hold on;
            subplot(1,2,2);
            [img, trks, colors]=highlightMq(img2, nz_mq, tail, k, Xh, c2, Xi, cmap);
            imshow(img); hold on;
        end

        for jj=nz_mq'
            r=getind(Xi.nX, k, jj, Xi.ri, 1);
            if mode==3
                trk=Xh(r,tail:k)/factor;
                plot3(trk(1,:), trk(2,:), trk(3,:), 'Color', cmap(jj,:), 'linewidth', 2);
                hold on;
                plot3(trk(1,end), trk(2,end), trk(3,end), 'o', 'Color', cmap(jj,:), 'linewidth', 2);
            else
                plot(trks(2*jj-1,:), trks(2*jj,:), 'Color', colors(jj,:));
                plot(trks(2*jj-1,end), trks(2*jj,end), 'Color', colors(jj,:));
    %             text(trks(2*jj-1,end)+3, trks(2*jj,end)+3, sprintf('%d', jj), 'Color', colors(jj,:));
            end
        end
        if mode==3, 
            axis(alim); 
            set(gca, 'fontsize', 16);
            box on;
            plotCompass(compass, alim, 16)
        end
        if ~record, title(k); end
        drawnow;
        if(record)
            set(gcf,'units','pixel');
    %         set(gcf,'position',[0,0,size(img1,2)/1.667,size(img1,1)/1.667]);
    %         set(gcf,'papersize',[size(img1,2)/1.667,size(img1,2)/1.667]);
            set(gcf,'PaperPositionMode','auto')
            %                 set(gca,'LooseInset',get(gca,'TightInset'));
            if mode==3
                set(gca, 'OuterPosition', [0 0 1 1]);
            else
    %             set(gca, 'OuterPosition', [-.18 -.18 1.35 1.35]);
            end
            print('-djpeg', sprintf('/Volumes/Puneet_Back/dumpss/stereo_%.4d.jpg', k), '-r200');
        end
        if ~record, pause(.1); end
    end
end

function [img, trks, colors]=highlightMq(img, nz_mq, tail, k, Xh, cam, Xi, cmap)

[h, w]=size(img(:,:,1));

% set the circular region for each mosquito
npts=200; th=linspace(0,2*pi, npts);
hl_rad=0:15;
[HLR, TH]=meshgrid(hl_rad, th);
regionx0=ceil(HLR.*cos(TH));
regiony0=ceil(HLR.*sin(TH));

trks=[];
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
    [num, thresh]=hist(double(img(imgind)),5);
    thresh=mean(thresh);
    img(imgind)=img(imgind)*1.1;
    for gg=1:numel(imgind)
        if img(imgind(gg))<thresh(1)
            img(imgind(gg))=img(imgind(gg))*.99;
        end
    end
end

function plotCompass(compass, als, fs)

    
compass=[compass(1,:)+((als(1,2)-als(1,1))/2);
         compass(2,:)+((als(1,4)-als(1,3))/2);  ];
% 
[val idx1]=min(compass(1,:));
[val idx2]=min(compass(2,:));
arrowtails=[compass(:,idx1), compass(:,idx2)];
% 
[val idx1]=max(compass(1,:));
[val idx2]=max(compass(2,:));
arrowheads=[compass(:,idx1), compass(:,idx2)]+.1;
% 
% % ------------- compass
% patch(compass(1,:), compass(2,:), ...
% als(1,5)*ones(1,size(compass,2)), ones(1,size(compass,2)), 'facecolor', ones(1,3)*.85);
quiver3(arrowtails(1,1), arrowtails(2,1), als(1,5), 1, 0, 0, .4, 'k');
% % quiver3(arrowtails(1,2), arrowtails(2,2), als(1,5), 0, 1, 0, .4, 'k');
% % text(arrowtails(1,1), arrowtails(2,1), als(1,5), 'S', 'fontsize', fs);
% % text(arrowtails(1,2), arrowtails(2,2), als(1,5), 'E', 'fontsize', fs);
text(arrowheads(1,1), arrowheads(2,1), als(1,5), 'N', 'fontsize', fs);
% text(arrowheads(1,2), arrowheads(2,2), als(1,5), 'W', 'fontsize', fs);
