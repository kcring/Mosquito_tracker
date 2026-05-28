function savename=cam2world(varargin)
% to convert the posfile into mat file for use in visualization and
% analysis

run config

if nargin==1
    % dataloc
%     pathname=varargin{1};
    % data_mq_XXX.mat
    filename=varargin{1};
    [pathname, name]=fileparts(filename);
    
    load(filename, 'Xh', 'Xi', 'P');
end
if ~strcmp(name(end-1:end), '_F') 
    fprintf('[!] Accepts smooth data only...\n');
    return;
end


savename=[pathname, '/',name, '_W.mat'];

% update floc anyway...
a=textscan(pathname, '%s', 'Delimiter', '/');
bb=a{:};
floc=sprintf('%s/', bb{1:end-2});

%% Processing
% Calibration information
addpath([floc '/calib']);
calibfun=dir([floc 'calib/' '*calib*.m']);
if(size(calibfun,1))
    get_cam_calib=str2func(calibfun.name(1:end-2));
else
    fprintf('[?] Could not find calibration for the cameras. Is this 2D tracking? ...\n');
end


%% convert to world coordinates
fprintf('Converting to world-frame coordinates....\n');
fprintf('Loading metadata from %s ....\n', floc);


metafile=dir([floc 'calib/' 'metadata_*.txt']);
if(size(metafile,1))
    meta_file=[floc, 'calib/', metafile(1).name];
else
    fprintf('[!]Could not find a metadata file\n');
    return;
end

meta_info=csvread(meta_file);
compass_dir=meta_info(1);
th_i = meta_info(2);
if numel(meta_info)<3
    fprintf('[!] Add camera height in mm to metadata in the 3rd column..\n');
    return
else
    fprintf('[I] Using camera height = %d mm \n', meta_info(3));
end
camera_height=meta_info(3);

% stop removing any zero rows... to maintain ids
% Xh=Xh(any(Xh,2),:);

% find nz rows
nz_rows=find(sum(Xh,2)~=0);
nz_row1=nz_rows(1);

nz_k_ind=find(Xh(nz_row1,:)~=0);
k0=nz_k_ind(1);

% swarm centroid
nt_k_ind=find(sum(Xh(1:Xi.nX:end,:),2)~=0);
% rc=[mean(Xh((nt_k_ind-1)*Xi.nX+1, k0)); mean(Xh((nt_k_ind-1)*Xi.nX+2,k0)); mean(Xh((nt_k_ind-1)*Xi.nX+3,k0))];
% T1=[eye(3), -rc;
%         0 0 0 1];
cam1=get_cam_calib(1);
cam2=get_cam_calib(2);
stereo_bar_center=(cam1.trm(1:3,4)+cam2.trm(1:3,4))/2;
% T1 moves to the center of the stereo bar
T1=[eye(3), -stereo_bar_center;
    0   0   0   1];

% T2 makes the world frame horizontal
th=pi/2-th_i*pi/180;   
T2=[inv(rota(th, 'x')), [0 0 0]';
       0 0 0 1];
   
% T3 rotates the world frame about the vertical (z) axis so that the x-axis is pointing 
% north (or y is pointing west)   
% compass dir is where the cameras are pointing. We want our z axis to
% always point towards the west 
th = (compass_dir-270)/180*pi;
T3=[inv(rota(th, 'z')), [0 0 0]';
       0 0 0 1];

% T4 drops the world frame to the ground right under the stereo bar center   
% T4=[eye(3), [0 0 1250]';
%         0 0 0 1];
T4=[eye(3), [0 0 camera_height]';
        0 0 0 1];


w1Tw0=T4*T3*T2*T1;
% to transform velocities in the inertial frame
w1Rw0=[ w1Tw0(1:3,1:3), [0 0 0]'
        0,  0,  0,  1];

%     nt=size(Xh,1)/Xi.nX;
for mq_ii = nt_k_ind'
    [r c]=getind(Xi.nX, 1, mq_ii, 1:Xi.nX, 1);
    % convert 
    trk=Xh(r(Xi.ri),:);
    nz_ind=find(trk(1,:)~=0);
    trk(:,nz_ind)=tra2b(trk(:,nz_ind), w1Tw0);
    Xh(r(Xi.ri),:)=trk;

    trk_vel=Xh(r(Xi.rdi),:);
    nz_ind=find(trk_vel(1,:)~=0);

    trk_vel(:,nz_ind)=tra2b(trk_vel(:,nz_ind), w1Rw0);
    Xh(r(Xi.rdi),:)=trk_vel;
end

save(savename, 'Xh', 'Xi');
% write csv files also as per jose's cycle analysis format
mq=1;
lastrow=0;
for jj=1:size(Xh,1)/Xi.nX
    r=getind(Xi.nX, 1, jj, 1:Xi.nX, 1);
    nzi=find(Xh(r(1),:)~=0);
    nframes=numel(nzi);
    if nframes
        cycle_analysis(lastrow+1:lastrow+nframes,:)=[nzi', jj*ones(nframes,1), Xh(r(1:6),nzi)'];
        lastrow=nframes+lastrow;
        mq=mq+1;
    end
end
csvwrite(sprintf('%s/%s_cycle_analysis_W.csv', pathname, name), cycle_analysis);


%for mq_ii= nt_k_ind'
%    [r c]=getind(Xi.nX, 1, mq_ii, 1:Xi.nX, 1);
%    if ~isempty(find(Xh(r(Xi.fi(2)),:)==2)) % only hand tracked mosquitoes
%        fprintf('Saving mq-%.2d_F_W.csv ....\n', mq_ii);
%        csvwrite(sprintf('%s/mq-%.2d_F_W.csv', pathname, mq_ii), Xh(r,:));
%    end
%end

% also save the new world coordinates
save([pathname, '/frame_W.mat'], 'w1Tw0');
