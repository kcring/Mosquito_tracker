function plot3Dfig

close all

if nargin < 1
    [floc, frmloc]=getexpdirs;
end

[optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams]=initproc(floc, frmloc);

[filename, pathname]=uigetfile('*.mat; *.csv', 'Pick a .csv file or a .mat file', ...
        [floc, 'output/data']);
[pathstr, name, ext]=fileparts(filename);
if strcmp(ext, '.mat')
    load([pathname, filename], 'Xh');
elseif strcmp(ext, '.csv')
    Xh=csvread([pathname, filename]);
end

if strcmp(pathname(end-11:end), 'output/data/')
    floc=pathname(1:end-12);
end

% find nonzero time-steps
Xh=Xh(any(Xh,2),:);
nz_k=find(Xh(1,:)~=0);


k0=nz_k(1);
kF=nz_k(end);
tl=5;
fs=24;

% convert everything to m
Xh=Xh/1000;

pvals_x=Xh(1:Xi.nX:end,:);
pvals_y=Xh(2:Xi.nX:end,:);
pvals_z=Xh(3:Xi.nX:end,:);
cmap=colormap(lines(size(Xh,1)/Xi.nX));
lw=1;

border=.25; %m

axis_limits=[min(pvals_x(pvals_x(:)~=0)),  min(pvals_y(pvals_y(:)~=0)), min(pvals_z(pvals_z(:)~=0)); 
             max(pvals_x(pvals_x(:)~=0)),  max(pvals_y(pvals_y(:)~=0)), max(pvals_z(pvals_z(:)~=0))];
axis_limits=axis_limits+[-border*ones(1,3); border*ones(1,3)];         

axis_limits=reshape(axis_limits, [1,6]);

nt_k_ind=size(Xh,1)/Xi.nX;
    
for t_ii=1:nt_k_ind
    [r c]=getind(Xi.nX, k0, t_ii, Xi.ri, 1);
    tXh=Xh(r,k0:kF);
    tXh(tXh==0)=nan;

    % smoothing...
    %for ii=Xi.ri
    %    tXh(ii,:)=sma(tXh(ii,:), 3);
    %end

    plot3(tXh(1,:), tXh(2,:), tXh(3,:), 'Color', cmap(t_ii,:), 'LineWidth', lw);
    hold on;
    plot3(tXh(1,end), tXh(2,end), tXh(3,end), '.', 'Color', cmap(t_ii,:), 'MarkerSize', 7)
end


% draw cameras
if 0
    load([floc, 'output/data/frame_W.mat']);
    w1Tw0(1:3,4)=w1Tw0(1:3,4)/1000;
    for cc=1:2
        cam=get_cam_calib(cc);
        cTw=cam.trm;
        cTw(1:3,4)=cTw(1:3,4)/1000;
        fov=[36, 32];
        
        drawCam(fov, 1, cTw*inv(w1Tw0), 1);
    end
else
    axis(axis_limits);
end


box on;
set(gca, 'FontSize', fs);
xlabel('m'); ylabel('m'); zlabel('m');
% quiver3(axis_limits(1), axis_limits(3)-.25, axis_limits(5), 1, 0, 0, 'k');
% text(axis_limits(1)+1, axis_limits(3)-.25, axis_limits(5), 'N', 'FontSize', 12);
drawnow;


print('-dpng', sprintf('%s/output/movies/fig_3dplots_%.4d-%.4d.png', floc, k0, kF));

pause(.05);


