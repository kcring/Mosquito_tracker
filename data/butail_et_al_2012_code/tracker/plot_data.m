function plot_data(datafile, frmloc, prntloc, k0, kF, mode)
%
% if only datafile is passed then 3D is plotted
% if frmloc is also passed then the estimates are reprojected
% on the screen
% 
% mode=1 plot all data on a 3D plot
% mode=2 write stereo images
if nargin < 1
    datafile = '~/Dropbox/EASeL/marker_exp/data/20160326_191454/output/data/data_mq_EXP_2016_03_26_L_19.14.mat';
    frmloc = '/Volumes/DUMP01/2016/nimr/20160326_191454/frames/';
    prntloc= frmloc;
    k0 = 488; kF = 546;
    mode=2; 
end
% 
% datafile = '/Users/puneet2895/Dropbox/marker_exp/data/20160326_190755/output/data/data_mq_EXP_2016_03_26_L_19.07.mat';
% frmloc = '/Volumes/Puneet_Back/20160326_190755/frames/';
% k0 = 385; kF = 585;
% 
% 
% datafile = '/Users/puneet2895/Dropbox/marker_exp/data/20160326_192011/output/data/data_mq_EXP_2016_03_26_L_19.14.mat';
% frmloc = '/Volumes/Puneet_Back/20160326_192011/frames/';
% k0 = 300; kF = 420;

% datafile = '/home/puneetjain/Dropbox/marker_exp/data/20160326_191707/output/data/data_mq_EXP_2016_03_26_L_19.17.mat';
% frmloc = '/media/puneetjain/Puneet_Back/20160326_191707/frames/';
% k0 = 187; kF = 687;

% datafile = '/Users/puneet2895/Dropbox/marker_exp/data/20160326_190755/output/data/data_mq_EXP_2016_03_26_L_19.07.mat';
% frmloc = '/Volumes/Puneet_Back/20160326_190755/frames/';
% k0 = 385; kF = 585;

% datafile = '/Users/puneet2895/Dropbox/marker_exp/data/20160325_191717/output/data/data_mq_EXP_2016_03_25_L_19.17.mat';
% frmloc = '/Volumes/Puneet_Back/20160325_191717/frames/';
% k0 = 15; kF = 165;

if mode==2
    load(datafile, 'Xh', 'Xi', 'cams', 'frmlist');
else
    load(datafile, 'Xh', 'Xi');
end

record=1;
% if record, visible=0; else visible=1; end

% ---- non-zero time-steps
nz=find(sum(Xh)~=0);
% k0=nz(1); kF=nz(end);
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
%         if visible
%             figure(1); gcf; clf;
%         else
%             figure('Visible', 'off'); clf;
%         end
        if mode==2
            img1=imread([frmloc, frmlist(k,1).name]); 
            img2=imread([frmloc, frmlist(k,2).name]); 
        end
        tail=max(k0, k-5);
        nz_mq=find(~isnan(Xh(1:Xi.nX:end,k)));


        if mode==2
            subplot(1,2,1);
            [img, trks1, colors1]=highlightMq(img1, nz_mq, tail, k, Xh, cams(1), Xi, cmap);
            imshow(img); hold on;
            subplot(1,2,2);
            [img, trks2, colors2]=highlightMq(img2, nz_mq, tail, k, Xh, cams(2), Xi, cmap);
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
                subplot(1,2,1);
                plot(trks1(2*jj-1,:), trks1(2*jj,:), 'Color', colors1(jj,:));
                plot(trks1(2*jj-1,end), trks1(2*jj,end), 'Color', colors1(jj,:));
                
                subplot(1,2,2);
                plot(trks2(2*jj-1,:), trks2(2*jj,:), 'Color', colors1(jj,:));
                plot(trks2(2*jj-1,end), trks2(2*jj,end), 'Color', colors1(jj,:));
    %             text(trks(2*jj-1,end)+3, trks(2*jj,end)+3, sprintf('%d', jj), 'Color', colors(jj,:));
            end
        end
        if mode==3 
            axis(alim); 
            set(gca, 'fontsize', 16);
            box on;
            plotCompass(compass, alim, 16)
        end
        if ~record, title(k); end
        drawnow;
        if(record)
            set(gcf,'units','normalized');
    %         set(gcf,'position',[0,0,size(img1,2)/1.667,size(img1,1)/1.667]);
    %         set(gcf,'papersize',[size(img1,2)/1.667,size(img1,2)/1.667]);
%             set(gcf,'PaperPositionMode','auto')
            set(gcf, 'position', [0.1677    0.7361    0.4589    0.1685]);
            subplot(1,2,1);
            s1p=get(gca, 'position');
            set(gca, 'position', [s1p(1)*1.8, s1p(2), s1p(3), s1p(4)]);
            subplot(1,2,2);
            s2p=get(gca, 'position');
            set(gca, 'position', [s2p(1)*0.9, s2p(2), s2p(3), s1p(4)]);
%             set(gca,'LooseInset',get(gca,'TightInset'));
%             if mode==3
%                 set(gca, 'OuterPosition', [0 0 1 1]);
%             else
%                 set(gca, 'OuterPosition', [-.18 -.18 1.35 1.35]);
%             end
            print('-djpeg', sprintf('%s/stereo_%.4d.jpg', prntloc, k), '-r250');
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
