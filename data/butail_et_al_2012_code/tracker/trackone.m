function varargout = trackone(varargin)
% TRACKONE M-file for trackone.fig
%      TRACKONE, by itself, creates a new TRACKONE or raises the existing
%      singleton*.
%
%      H = TRACKONE returns the handle to a new TRACKONE or the handle to
%      the existing singleton*.
%
%      TRACKONE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKONE.M with the given input arguments.
%
%      TRACKONE('Property','Value',...) creates a new TRACKONE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mqtrack_manual_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackone_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackone

% Last Modified by GUIDE v2.5 17-Jun-2012 12:18:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackone_OpeningFcn, ...
                   'gui_OutputFcn',  @trackone_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before trackone is made visible.
function trackone_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackone (see VARARGIN)

% --------- Set proper environment ------------------
if size(varargin,2) ~=2
    [floc, frmloc]=getexpdirs;
else
    floc=varargin{1};
    frmloc=varargin{2};
end

[optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams]=initproc(floc, frmloc);

% Trackone options
optTrack.trackone.alg='mpf'; % algorithm is 'ukf', 'pf', 'ekf'
optTrack.trackone.v0=[10 10 10]'; % mm/s
optTrack.trackone.dn=50; % mm (disturbance noise) standard deviation
optTrack.trackone.N=1000; % number of samples for particle filter
optTrack.trackone.mm_pix_t=10; % pixel radius while looking for temporal da
optTrack.trackone.camstd=[4 4]; % noise covariance for each cam
optTrack.trackone.threed_t=25; % mm distance where to look for auto-tracks.
optTrack.trackone.max_t=100;
optTrack.trackone.swarm_d0=swarm_boundaries(1); % mm
optTrack.trackone.swarm_s0=swarm_boundaries(2); % mm

% ---------- Pass variables ----------------------
handles.floc=floc;
handles.frmloc=frmloc;
handles.dataloc=dataloc;
handles.expname=expname;
handles.k0=1;
handles.k=handles.k0+optTrack.br0;
handles.get_cam_calib=get_cam_calib;
handles.cams=readOffCamCalib(calibloc, nc);
handles.optTrack = optTrack;
handles.searchbbox=optTrack.bbox;

% ------------- Image streams / camera specific info -------------
% camera identifiers (change as per file-naming conventions)
camids=['L'; 'R'];
handles.imageIds=imageIds;
handles.nc=nc;
% Get list of all image files
fprintf('Listing files ...\n');
for ii=1:handles.nc
    handles.imstream(ii).flist = dir([handles.frmloc, ...
                                handles.imageIds(ii,:), '*.*']);
    if(~size(handles.imstream(ii).flist,1))
        error('[!] Image files not found');
    end
end
% initialize thresholds
if ~exist(sprintf('%scalib/cam1_bgparams.csv', floc), 'file')
    expected_nmq=input('Max # of mosquitoes you expect to see in this swarm (approx.): ');
else
    expected_nmq=40;
end
for ii=1:size(camids,1)
    handles.imstream(ii).img_info=imfinfo([handles.frmloc, ...
                    handles.imstream(ii).flist(1).name]);
    %******MATLAB doesn't read the bitdepth properly ??%
    if handles.imstream(ii).img_info.BitDepth >16, handles.imstream(ii).img_info.BitDepth=8; end     
    handles.imstream(ii).bitval=2^handles.imstream(ii).img_info.BitDepth-1;               
    handles.imstream(ii).imgarr=zeros(handles.imstream(ii).img_info.Height, ...
                                handles.imstream(ii).img_info.Width, ...
                                handles.optTrack.br0*2+1);   
    bgparams=init_bgparams(handles.imstream(ii), ii, floc, frmloc, ...
                    expected_nmq, optTrack);
    handles.imstream(ii).binary_t=bgparams.binary_t;
    handles.imstream(ii).area_t=bgparams.area_t;
    handles.imstream(ii).br=bgparams.br;
    handles.imstream(ii).noise_std=bgparams.noise_std;                            
    handles.imstream(ii).roi=bgparams.roi;
end


% --------------- Load data files -------------------
nframes=min(size(handles.imstream(1).flist,1),size(handles.imstream(2).flist,1));
handles.datfile0 = [handles.dataloc, sprintf('data_mq_auto_%s.mat', expname)];
if exist(handles.datfile0, 'file')
    fprintf('[I] Found automatically tracked data... loading...\n');
    reformat_data(handles.datfile0);
    load(handles.datfile0, 'Xh', 'Xi', 'frmlist');
    m_mqid_all_auto=show_tracked_mq(Xh, Xi, 0, 0);
    mqmax=max(m_mqid_all_auto);
    handles.Xi0=Xi;
    handles.Xh0=Xh(1:mqmax*Xi.nX,:);
    handles.frmlist=frmlist;
    if exist('frames', 'var')
        handles.frame=frames;
    end
    clear('Xh');
else
    fprintf('[A] Proceeding without any auto tracks...\n');
    Xi=strXi;
    handles.Xi=Xi;
    handles.frmlist=[];
end

handles.p=zeros(Xi.nX, optTrack.trackone.N); % for particle filter
handles.datfile1 = [handles.dataloc, sprintf('data_mq_%s.mat', expname)];

if(~exist(handles.datfile1, 'file'))
    fprintf('Creating %s ... \n', handles.datfile1);
    Xh=zeros(Xi.nX*optTrack.trackone.max_t, nframes);
%     Z=strZ([optTrack.trackone.max_t,nframes]);
%     save(handles.datfile1, 'Xh', 'Z', 'Xi', 'optTrack', 'cams');
    save(handles.datfile1, 'Xh', 'Xi','cams');
end
reformat_data(handles.datfile1);
load(handles.datfile1, 'Xh');
handles.Xh=Xh;
% if ~isstruct(Z)
%     error('[!] This dataset was generated using the older version... run xx to reformat it\n');
% end
% handles.Z=Z;  
% this Xi value takes precendence over the earlier one
handles.Xi=Xi;

%%%% this part adds rows if not present
if size(handles.Xh, 1)/handles.Xi.nX < optTrack.trackone.max_t
    addz=(optTrack.trackone.max_t-size(handles.Xh, 1)/handles.Xi.nX)*Xi.nX;
    handles.Xh=[handles.Xh; ones(addz,1)*handles.Xh(1,:)*0];
end

%%%% this part adds columns if needed
if size(handles.Xh, 2) < nframes
    addc=nframes-size(handles.Xh, 2);
    handles.Xh=[handles.Xh, (handles.Xh(:,1)*0)*ones(1,addc)];
end


% ------------- Initialize other variables ---------------

handles.backtrack=0;
handles.mqid=0;

% GUI Look and Feel
handles.ms=[7 7]; % marker size on each frame
handles.lw=[1 1]; % line width on each frame
handles.cmap=colormap(jet(99));
handles.fz =1; % focused zoom factor
% show data
% find the mosquitoes that are being tracked manually at this time-step
m_mqid=show_tracked_mq(handles.Xh, handles.Xi, 0, 1);

[handles.mqid handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);
logger([handles.floc, '/output/trkrun.log'], sprintf('trackone.m, k=%d, mqid selected=%d', ...
        handles.k, handles.mqid));

% Show input data
fprintf('Data location = "%s\n"', handles.dataloc);
fprintf('Image identifier = "%s"\n', image_id)
fprintf('Calibration function = "%s"\n', func2str(handles.get_cam_calib));
fprintf('Loading existing data .....\n');


clear('Xh');


for ii=1:handles.nc
    handles=store_orig_image(ii, handles);
    
    % Call function to segment image
    handles=init_img_array(handles, ii);                
    handles=apply_bg_subtract(handles, ii);
end


fprintf('Total tracks=%d ...\n', numel(m_mqid));

logger([handles.floc, '/output/trkrun.log'], sprintf('trackone.m, start, floc=%s, image_id=%s', handles.floc, image_id));

% setting the frame number
set(handles.pb_goto_frm, 'String', sprintf('Goto # [%d - %d]', optTrack.br0+1, size(handles.Xh,2)-optTrack.br0));

% axes handles
handles.ah(1)=handles.axes1;
handles.ah(2)=handles.axes2;
handles.ah(3)=handles.axes3;
handles.ah(4)=handles.axes6;

% Default for magnification
for ii=1:handles.nc
    handles.imstream(ii).xlim=[0 handles.imstream(ii).img_info.Width];
    handles.imstream(ii).ylim=[0 handles.imstream(ii).img_info.Height];
end

% Set the L-R radio button on top of tuners to L
set(handles.rb_tune_R, 'Value', 1);
set(handles.rb_original_R, 'Value', 1);
set(handles.rb_original_L, 'Value', 1);
set(handles.cb_see_swarm, 'Value', 0);

% which camera has noise/area tuners by default
handles.tunecam=2;

set(handles.ed_frame, 'String', sprintf('%d', handles.k));

% Display image frame in the axis for each camera
for ii=1:handles.nc
    show_frames(ii, 0, handles);
end

% plot the speed and position
plot_curr_mq_k_3d_speed(handles);


% Remove the axes labels from each plot
for ii=1:handles.nc
    axes(handles.ah(ii));
    axis off;
end

% Choose default command line output for trackone
handles.output = hObject;

% keyboard navigation
set(handles.figure1, 'KeyPressFcn',@keypress_nav);

% Update handles structure
guidata(hObject, handles);

% set_instr('[I] Click on a mosquito on the left frame', handles)

% UIWAIT makes trackone wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function set_instr(str, handles)
% fprintf('[%d]%s\n', handles.k, str)

if strcmp(str(1:3), '[!]')
    set(handles.txt_instr1, 'BackgroundColor', 'r');
elseif strcmp(str(1:3), '[A]')
    set(handles.txt_instr1, 'BackgroundColor', 'b');
else
    set(handles.txt_instr1, 'BackgroundColor', [0 .75 0]);
end

strold=get(handles.txt_instr1, 'String');

str=sprintf('%s [%d, %d] %s;', strold, handles.k, handles.mqid, str);

str_cell=textscan(str, '%s', 'Delimiter', ';');

set(handles.txt_instr1, 'String', sprintf('%s;', str_cell{1}{end-1:end}));


% Optimization function                        
function f=undist(xyn, xyd, kc1, kc2)
rn2=xyn(1)^2+xyn(2)^2;
f=norm(xyd-xyn*(1+kc1*rn2 +kc2*rn2^2));

% --- Outputs from this function are returned to the command line.
function varargout = trackone_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sl_bgt_Callback(hObject, eventdata, handles)
% hObject    handle to sl_bgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% find max value for each pixel through all images
camid=handles.tunecam;
handles.imstream(camid).binary_t=get(hObject,'Value');
set(handles.txt_bin_t, 'String', sprintf('Intensity:%.4f', handles.imstream(camid).binary_t));

set(gcf,'Pointer','watch');
handles=extract_pts(handles, camid);

show_frames(camid, get_frame_choice(handles,camid), handles);
% if(camid==1 && get(handles.rb_segmented_L, 'Value'))
%     show_frames(camid, 1, handles);
% elseif(camid==2 && get(handles.rb_segmented_R, 'Value'))
%     show_frames(camid, 1, handles);
% end

% Update handles structure
guidata(hObject, handles);

set(gcf,'Pointer','arrow');
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



% --- Executes during object creation, after setting all properties.
function sl_bgt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_bgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in pb_prev1.
function pb_prev1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_prev1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pb_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.k-1) <= handles.optTrack.br0
    set_instr('[!] You are on the earliest possible frame', handles)
else

    handles.k = handles.k - 1;

    set(handles.ed_frame, 'String', sprintf('%d', handles.k));


    handles.backtrack=1;

    % Store all images that are needed to segment
    for ii=1:handles.nc
        handles=store_orig_image(ii, handles);
        % Call function to segment image
        if(get(handles.cb_track, 'Value'))
            handles=update_img_array(handles, ii);
            handles=apply_bg_subtract(handles, ii);
        end
    end


    % Get axes values for magnification
    for ii=1:handles.nc
        axes(handles.ah(ii));
        handles.imstream(ii).xlim=get(gca, 'xlim');
        handles.imstream(ii).ylim=get(gca, 'ylim');
    end

    % track backwards also ....
    [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1:handles.Xi.nX, 1);
    if(get(handles.cb_track, 'Value'))
        if(~handles.Xh(r(1),c)|| get(handles.cb_override_existing_track, 'Value'))
            handles = predict(handles);
        end
    end


    show_frames(1,get_frame_choice(handles,1), handles);


    if(get(handles.cb_track, 'Value'))
        if(~handles.Xh(r(1),c)|| get(handles.cb_override_existing_track, 'Value'))
            handles=update(handles);
%         else
%             set_instr('[A] Click on a mosquito in the left frame', handles)
        end
    end

    show_frames(2,get_frame_choice(handles,2), handles);

    for c_ii=1:handles.nc
        plot_curr_mq_k(handles, c_ii);
        if get(handles.cb_see_swarm, 'Value')
            plot_all_mq_k(handles, c_ii);
        end
    end

    % start change 1/24/2011 : show 3D position and speed
    plot_curr_mq_k_3d_speed(handles);
    % end change 1/24/2011 : show 3D position and speed
    
    % set_instr('Click on the mosquito (or possible choices) in the left frame', handles)
    guidata(hObject, handles);
end

% set_instr('Click on the mosquito (or possible choices) in the left frame', handles)

% guidata(hObject, handles);


function plot_all_mq_k(handles, camid)
% plot all points available in current frame
nz_ind=find(handles.Xh(1:handles.Xi.nX:end,handles.k)~=0);
nt=0;
for ii=nz_ind'
    nt=nt+1;
    [r, c]=getind(handles.Xi.nX, handles.k, ii, handles.Xi.ri, 1);
    pos_3d(:,nt)=handles.Xh(r, c);
end

if(~isempty(nz_ind))

    pix=w2cam(pos_3d, handles.get_cam_calib(camid));
    axes(handles.ah(camid));

    for t_ii=1:nt
        plot(pix(1,t_ii), pix(2,t_ii), 's', 'Color', handles.cmap(t_ii,:), 'MarkerSize', 3);
        text(pix(1,t_ii)+rand*3+1, pix(2,t_ii), sprintf('%d', nz_ind(t_ii)), 'Color', handles.cmap(t_ii,:));
    end
    
    for t_ii=1:nt
        [r, c]=getind(handles.Xi.nX, handles.k, ii, handles.Xi.ri, 1);
        tr_3d=handles.Xh(r, :)/1000;
        tr_3d(:,tr_3d(3,:)==0)=nan;
        tr_3d(:,1:handles.k-5)=nan;
        tr_3d(:,handles.k+3:end)=nan;
        plot3(handles.ah(3), tr_3d(1,:), tr_3d(2, :),tr_3d(3, :),...
            '-', 'Color', handles.cmap(t_ii,:));
    end
end
drawnow;

function plot_curr_candidates(handles, camid)
jj=handles.mqid;
axes(handles.ah(camid));
cc=camid;
cmap_id=find(handles.mqids==handles.mqid);
% ya ya... I know lot of ifs will fix it soon.... 1/11
if isfield(handles.frame(handles.k), 'mq')
    if size(handles.frame(handles.k).mq, 2) >=jj
        if isfield(handles.frame(handles.k).mq(jj), 'cam')
            if size(handles.frame(handles.k).mq(jj).cam, 2) >=cc
                if isfield(handles.frame(handles.k).mq(jj).cam(cc), 'candpos')
                    if ~isempty(handles.frame(handles.k).mq(jj).cam(cc).candpos)
                        [jk, nt]=size(handles.frame(handles.k).mq(jj).cam(cc).candpos);
                        for ii=1:nt
                            plot(handles.frame(handles.k).mq(jj).cam(cc).candpos(1,ii), ...
                                                 handles.frame(handles.k).mq(jj).cam(cc).candpos(2,ii), ...
                                                 'x', 'Color', handles.cmap(cmap_id,:)+(nt-ii)/nt*([1 1 1]-handles.cmap(cmap_id,:)), ...
                                                 'MarkerSize', handles.ms(cc));
                                             hold on;
                        end
                    end
                end
            end
        end
    end
end

% start change 1/24/2011 : show 3D position and speed
function plot_curr_mq_k_3d_speed(handles)

[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1:handles.Xi.nX, 1);
tl=5;

t0=handles.k-tl;
if t0 < 1, t0=1; end

tf=handles.k+tl;
if tf > size(handles.imstream(1).flist,1)-handles.imstream(1).br 
    tf=size(handles.imstream(1).flist,1)-handles.imstream(1).br; 
end

% check if there is any value in the current timestep

if(size(handles.Xh,1)>=r(1))
    if sum(handles.Xh(r(1),:))
        tXh=handles.Xh(r, t0:tf);
        tXh(tXh==0)=nan;

        cmap_id=find(handles.mqids==handles.mqid);

        % ------ 3D Plot --------------
        axes(handles.ah(3));
        cla;
        pos=tXh(handles.Xi.ri,:);
        plot3(pos(1,:)/1000, pos(2,:)/1000, pos(3,:)/1000, '-', 'Color', handles.cmap(cmap_id,:));
        hold on;
        plot3(pos(1,:)/1000, pos(2,:)/1000, pos(3,:)/1000, '.', 'Color', handles.cmap(cmap_id,:));    
        plot3(pos(1,tl+1)/1000, pos(2,tl+1)/1000, pos(3,tl+1)/1000, 'o', 'Color', ...
                    handles.cmap(cmap_id,:), 'MarkerSize', 7);
        grid on;

        xlabel('x(m)'); ylabel(handles.ah(3), 'y(m)');
        set(gca, 'FontSize', 9);

        % --------- speed ----------
        axes(handles.ah(4));
        cla;
        vel=tXh(handles.Xi.rdi,:);
        speed=sqrt(sum(vel.^2));
        plot(t0:tf, speed/1000, 'Color', handles.cmap(cmap_id,:));
        grid on;
        hold on;
        plot(handles.k, speed(tl+1)/1000, 'o', 'Color', handles.cmap(cmap_id,:), 'MarkerSize', 7);
        % warning line.. at 4 m/s
        plot([t0 tf], [3.5 3.5], 'r', 'LineWidth', 2);
        set(gca, 'Ylim', [0, 4.5]);
        set(gca, 'Xlim', [t0 tf]);
        ylabel('speed (m/s)');
        xlabel('Frame #');

        set(gca, 'FontSize', 9);
    end    
end

drawnow;
% end change 1/24/2011 : show 3D position and speed

function plot_curr_mq_k(handles, camid)

r=getind(handles.Xi.nX, handles.k, handles.mqid, handles.Xi.ri, 1);
tl=5;

t0=handles.k-tl;
if t0 < 1, t0=1; end

tf=handles.k+tl;
if tf > size(handles.imstream(1).flist,1)-handles.imstream(1).br 
    tf=size(handles.imstream(1).flist,1)-handles.imstream(1).br; 
end

% check if there is any value in the current timestep
% if the number of rows
if(size(handles.Xh,1)>=r(1))
    if sum(handles.Xh(r(1),:))
        tXh=handles.Xh(r, t0:tf);
        tXh(tXh==0)=nan;


        pix=w2cam(tXh, handles.get_cam_calib(camid));
        cmap_id=find(handles.mqids==handles.mqid);


%         axes(handles.ah(camid));
        if(get(handles.cb_trajectory, 'Value'))        
            plot(handles.ah(camid), pix(1,:), pix(2,:),  '-o', 'Color', ...
                        handles.cmap(cmap_id,:)+.5*([1 1 1]-handles.cmap(cmap_id,:)), 'MarkerSize', handles.ms(camid));
            if(handles.backtrack)
                plot(handles.ah(camid), pix(1,1), pix(2,1), '.', 'Color', [.5 .5 .5], 'MarkerSize', handles.ms(camid));        
            else
                plot(handles.ah(camid), pix(1,end), pix(2,end), '.', 'Color', [.5 .5 .5], 'MarkerSize', handles.ms(camid));                
            end

            plot(handles.ah(camid), pix(1,tl+1), pix(2,tl+1), 'o', 'Color', handles.cmap(cmap_id,:)+.25*([1 1 1]-handles.cmap(cmap_id,:)),...
                    'MarkerSize', handles.ms(camid)*2);
        else
            % make it darker if you are not watching the trajectory
            plot(handles.ah(camid), pix(1,tl+1), pix(2,tl+1), '.', 'Color', handles.cmap(cmap_id,:),...
                'MarkerSize', handles.ms(camid));
        end
    end
end

drawnow;
    
function handles = predict(handles)
pr_mot=[];
jj=handles.mqid;

% get lcam structure
lcam=handles.get_cam_calib(1);

% last predicted target position
if(handles.backtrack)
    kt=handles.k+1;
else
    kt=handles.k-1;
end
[r, c]=getind(handles.Xi.nX, kt, handles.mqid, 1:handles.Xi.nX, 1);
rprev=handles.Xh(r(handles.Xi.ri),c);
vprev=handles.optTrack.trackone.v0;

% update velocity estimate if available
if(size(handles.Xh,1)>=r(handles.Xi.rdi(1)))
    if(handles.Xh(r(handles.Xi.rdi(1)),c))
        vprev=handles.Xh(r(handles.Xi.rdi),c);
    end
end

% ** get current object. if it was the prev button then make vprev negative
if (handles.backtrack)
    vprev=-vprev;
end

switch handles.optTrack.trackone.alg
    case 'mpf'
    % sample from a 3 dimensional gaussian pdf
    wk=randn(handles.optTrack.trackone.N,3)*handles.optTrack.trackone.dn;
    wk_pdf=normpdf(wk(:,1), 0, handles.optTrack.trackone.dn).*...
        normpdf(wk(:,2), 0, handles.optTrack.trackone.dn).*...
        normpdf(wk(:,3), 0, handles.optTrack.trackone.dn);

    % motion model (N particles in the direction of motion)
    next_pos_3d=repmat(rprev+handles.optTrack.dt*vprev, 1, handles.optTrack.trackone.N)+...
        (eye(3)*handles.optTrack.dt^2/2*wk');
    rkm=mean(next_pos_3d,2);
    lcam_next_pos_2d=w2cam(next_pos_3d, lcam);
    lcam_next_pos_2d_mean=w2cam(rkm, lcam);
    
    % find the pdf sum for each of these points and a point on the left
    % frame; a circle where you think the points may lie will help reduce
    % the number of points to evaluate

    lcam_dots=handles.imstream(1).centroids;
    nc1=size(lcam_dots,1);

    zz=0;
    for ii=1:nc1
        % if the centroid in the left image is within a threshold distance
        if(sqrt((lcam_dots(ii,1)-lcam_next_pos_2d_mean(1))^2+...
                (lcam_dots(ii,2)-lcam_next_pos_2d_mean(2))^2)<handles.optTrack.trackone.mm_pix_t)
            zz=zz+1;
            % likelihood of each of the dots w.r.t projected dots
            lcam_pdf=normpdf(lcam_next_pos_2d(1,:), lcam_dots(ii,1), handles.optTrack.trackone.camstd(1)).*...
                     normpdf(lcam_next_pos_2d(2,:), lcam_dots(ii,2), handles.optTrack.trackone.camstd(1));
            lcam_selected_dots(:,zz)=lcam_dots(ii,:)';
            % probability of the dots times the probability of the position
            % itself (from disturbance)
            pr_mot(zz)=sum(lcam_pdf.*wk_pdf');
        end
    end

    axes(handles.ah(1))
    if(~isempty(pr_mot))
        [val, idx]=sort(pr_mot);
        num_pr=length(pr_mot);

        %choose top 5 or num_pr whichever is less
        if(num_pr <5)
            handles.frame(handles.k).mq(jj).cam(1).candpos=lcam_selected_dots(:,idx);
        else
            handles.frame(handles.k).mq(jj).cam(1).candpos=lcam_selected_dots(:,idx(end-4:end));
        end

        fprintf('[%d, %d]\n',handles.k, jj);
        % get ready to click on a mosquito
        handles.frame(handles.k).mq(jj).c(:,1)=lcam_selected_dots(:,idx(end));
%         [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.uv, 2);
%         handles.Z(r,c(1))=lcam_selected_dots(:,idx(end));
        handles.Z(2*handles.mqid-1,handles.k)=handles.imstream(1).Z(idx(end));
    else
          set_instr('[!] Could not find any targets within the probable region... select manually', handles);
%           pause(2);
          handles.frame(handles.k).mq(jj).c(:,1)=lcam_next_pos_2d_mean;
%           [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.uv, 2);
%           handles.Z(r,c(1))=lcam_next_pos_2d_mean;
    end
    case 'pf'
        
        handles.p(r,:)=mq_motion(handles.p(r,:), optTrack);
end

% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ua, 2);
% handles.Z(r,c(1))=1;

%---------- commenting
% for ii=1:handles.nc
%     plot_curr_candidates(handles, ii);
% end



function dummy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_bgt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in cb_lc.
function cb_lc_Callback(hObject, eventdata, handles)
% hObject    handle to cb_lc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_lc

if(get(hObject,'Value'))
    handles=update(handles);
end


guidata(hObject, handles);


function handles=update(handles)

Pr=[];
jj=handles.mqid;

% Cylinder points 
% To be changed to a truncated cone later
rc=0:50;
tc=linspace(-pi, pi, 30);

startp=handles.optTrack.trackone.swarm_d0;
endp=handles.optTrack.trackone.swarm_d0+handles.optTrack.trackone.swarm_s0;

if(handles.backtrack)
    kt=handles.k+1;
else
    kt=handles.k-1;
end
[r, c]=getind(handles.Xi.nX, kt, handles.mqid, handles.Xi.ri, 1);

if(handles.Xh(r(1),c))
%     startp=norm(handles.frame(handles.k-1).mq(jj).X(1:3))-300;
    startp=norm(handles.Xh(r,c))-300;
    endp=startp+500;
end

switch handles.optTrack.trackone.alg
    case 'mpf'
        dc=0:5:(endp-startp);
        [Rc, Tc, Dc]=meshgrid(rc,tc,dc);

        Xc=Rc.*cos(Tc);
        Yc=Rc.*sin(Tc);
        Zc=Dc;

        %%% projecting epipolar line %%

        lcam=handles.get_cam_calib(1);
        rcam=handles.get_cam_calib(2);

        % store mosquito values in a local variable
        if(isfield(handles, 'frame'))
            if(size(handles.frame, 2)>=handles.k)
                lc=handles.frame(handles.k).mq(jj).c(:,1);
            end
        end

        axes(handles.ah(1));

        if(exist('lc', 'var'))
            eline_w=get_eline([lc(1), lc(2)], ...
                        lcam, startp, endp);


            eline_rcam=w2cam(eline_w, rcam);   
%             handles.frame(handles.k).reline(:,:,jj)=eline_rcam;
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
            c1=c0-handles.optTrack.bbox(1);
            c3=c0+handles.optTrack.bbox(1);

            % swarm size in pixels
            ssp=sqrt((eline_rcam(1,end)-eline_rcam(1,1))^2 + ...
                     (eline_rcam(2,end)-eline_rcam(2,1))^2);
            c2=cp;
            c4=-eline_rcam(1,end)*mp + eline_rcam(2,end);

            % number of dots in right cam
            rcam_dots=handles.imstream(2).centroids;

            nce=size(rcam_dots,1);
            eline_v3=(eline_w(:,2)-eline_w(:,1))/norm(eline_w(:,2)-eline_w(:,1));
            % making an orthogonal frame since the directions of other two don't
            % matter we choose v1=[a b c], v2=[0 -c -b] and v3= v1xv2;
            eline_v2=[0 -eline_v3(3) eline_v3(2)]'/norm([0 -eline_v3(3) eline_v3(2)]);
            eline_v1=cross(eline_v2, eline_v3);
            wTeline=[eline_v1, eline_v2, eline_v3, eline_w(:,1);
                        0 0 0 1];

            % points in 3D where we need to find probability values
            [re1, re2, re3]=trpa2b(Xc, Yc, Zc, wTeline);
            [uCL, vCL]=wp2cam(re1, re2, re3, lcam);
            Pr_camL=normpdf(uCL, lc(1), handles.optTrack.trackone.camstd(1)).*...
                    normpdf(vCL, lc(2), handles.optTrack.trackone.camstd(1));

            [uCR, vCR]=wp2cam(re1, re2, re3, rcam);

            zz=0;
            selected_dots_idx=[];
            for ii=1:nce
                xi=rcam_dots(ii,1);
                yi=rcam_dots(ii,2);

                % to fit in a tilted box we ask that the y values be between the
                % enclosing lines m1 is the line parallel to epipolar line just
                % below it (by handles.optTrack.bbox) m2 is normal to it at one end, m3 is
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
                    selected_dots_idx=[ii, selected_dots_idx];
                end
            end

            for ii=1:zz
                xi=rcam_selected_dots(1,ii);
                yi=rcam_selected_dots(2,ii);
                Pr_camR=normpdf(uCR, xi, handles.optTrack.trackone.camstd(2)).*...
                    normpdf(vCR, yi, handles.optTrack.trackone.camstd(2));
                pdf=Pr_camL.*Pr_camR;
                [valm, idxm]=max(pdf(:));
                [i1, i2, i3]= ind2sub(size(pdf), idxm);
                map_r(:,ii)=[re1(i1,i2,i3), re2(i1,i2,i3), re3(i1,i2,i3)]';
                Pr(ii)=sum(sum(sum(pdf)));    
            end

            handles.frame(handles.k).mq(jj).eline=eline_rcam;
            plot_eline(handles, 2);

            if(~isempty(Pr))
                [val, idx]=sort(Pr);
                num_Pr=length(Pr);

                % choose top 5 or num_Pr whichever is less
                if (num_Pr <5)
                    handles.frame(handles.k).mq(jj).cam(2).candpos=rcam_selected_dots(:,idx);
                else
                    handles.frame(handles.k).mq(jj).cam(2).candpos=rcam_selected_dots(:,idx(end-4:end));
                end

                    handles.frame(handles.k).mq(jj).c(:,2)=rcam_selected_dots(:,idx(end));

%                     [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.uv, 2);
%                     handles.Z(r,c(2))=rcam_selected_dots(:,idx(end));
%                     handles.Z(2*handles.mqid, handles.k)=handles.imstream(2).Z(selected_dots_idx(idx(end)));
%                     [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ua, 2);
%                     handles.Z(r,c(2))=1;

                    [r, c]=getind(handles.Xi.nX, handles.k, jj, 1:handles.Xi.nX, 1);
                    handles.Xh(r(handles.Xi.ri), c)=map_r(:,idx(end));

                    % *** if backtracking 
                    if (handles.backtrack)
                        if(handles.Xh(r(1), c+1))
                            handles.Xh(r(handles.Xi.rdi),c)=(handles.Xh(r(handles.Xi.ri),c+1)-handles.Xh(r(handles.Xi.ri),c))/handles.optTrack.dt;
                         end
                    else
                        if(handles.Xh(r(1), c-1))
                            handles.Xh(r(handles.Xi.rdi),c)=(handles.Xh(r(handles.Xi.ri),c)-handles.Xh(r(handles.Xi.ri),c-1))/handles.optTrack.dt;
                        end
                    end
                    % nf(2) == 2 for manual tracking
%                     handles.Xh(r(handles.Xi.fi), c)=[-1; 2];
                    instr=sprintf('3D position: [%d, %d, %d]\nPress next (>) or save', ...
                                    ceil(map_r(1,idx(end))), ceil(map_r(2,idx(end))), ceil(map_r(3,idx(end))));
                    set_instr(instr, handles)

            else
                set_instr('[!] No match found in the right camera image... Mark a point manually (also check binary threshold)', handles);
                % 
            end
        end
    case 'pf'
        
        
end

% call splice tracks after every update
% change 9/5/2018 Sachit - add override to the condition
if isfield(handles, 'Xh0') && ~get(handles.cb_override_existing_track, 'Value')
    [handles.Xh, ind]= t1_splice_tracks(handles.Xh, handles.Xh0, handles.mqid, handles.k, handles.backtrack,...
                    handles.Xi, handles.Xi0, handles.optTrack.trackone.threed_t);
    if ind
        set_instr('[A] Track spliced to an auto version... Check override or keep moving', handles);
    end                
end


% for ii=1:handles.nc
%     plot_curr_candidates(handles, ii);
% end

function plot_eline(handles, camid)

axes(handles.ah(camid));
eline_cam=handles.frame(handles.k).mq(handles.mqid).eline;
plot(eline_cam(1,:), eline_cam(2,:), 'Color', [.65 .65 .65]);


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

function handles = store_orig_image(camid, handles)

kk=handles.k;
imglist=handles.imstream(camid).flist;
try
    zz=imread([handles.frmloc, imglist(kk).name]);
    handles.imstream(camid).orig_img=zz;
catch
    error('Image file not found ... ');
end
if(camid==1)
    set(handles.txt_l_filename, 'String', imglist(kk).name);
elseif(camid==2)
    set(handles.txt_r_filename, 'String', imglist(kk).name);
end

function handles = init_img_array(handles, camid)
% call this function every time we jump to a new time-step
jj1=1;

for jj=handles.k-handles.optTrack.br0:handles.k+handles.optTrack.br0
    img0=imread([handles.frmloc, handles.imstream(camid).flist(jj).name]);
    if(size(img0,3)>1), img0=rgb2gray(img0); end
    if(handles.imstream(camid).noise_std>0), img0=filter2(fspecial('gaussian', [3 3], handles.imstream(camid).noise_std), img0); end
    % index is such that we want k to be right on optTrack.br0+1,
    % therefore we go back br0 and add jj and then -1
    handles.imstream(camid).imgarr(:,:,jj1)=img0;
    jj1=jj1+1;    
end

function handles = update_img_array(handles, camid)

if(handles.backtrack)
    handles.imstream(camid).imgarr=circshift(handles.imstream(camid).imgarr, [0,0, 1]);
    % populate the last index with k+optTrack.br0 image
    img0=imread([handles.frmloc, handles.imstream(camid).flist(handles.k-handles.optTrack.br0).name]);
    if(size(img0,3)>1), img0=rgb2gray(img0); end
    if(handles.imstream(camid).noise_std>0), img0=filter2(fspecial('gaussian', [3 3], handles.imstream(camid).noise_std), img0); end        
    handles.imstream(camid).imgarr(:,:,1)=img0;
else
    handles.imstream(camid).imgarr=circshift(handles.imstream(camid).imgarr, [0,0, -1]);
    % populate the last index with k+optTrack.br0 image
    try
        img0=imread([handles.frmloc, handles.imstream(camid).flist(handles.k+handles.optTrack.br0).name]);
    catch ME
        fprintf('[!] Could not find the frame....');
        throw(ME);
    end
    if(size(img0,3)>1), img0=rgb2gray(img0); end
    if(handles.imstream(camid).noise_std>0), img0=filter2(fspecial('gaussian', [3 3], handles.imstream(camid).noise_std), img0); end        
    handles.imstream(camid).imgarr(:,:,2*handles.optTrack.br0+1)=img0;
end

function handles = apply_bg_subtract(handles, camid)

% [binary_t area_t]=t1_retrieve_thresholds(handles.Z, handles.k, handles.mqid, camid, handles.Zi);
% 
% if binary_t, handles.imstream(camid).binary_t=binary_t; end
% if area_t, handles.imstream(camid).area_t=area_t; end


set(handles.sl_bgt, 'Value', handles.imstream(camid).binary_t);
set(handles.txt_bin_t, 'String', sprintf('Intensity:%.4f', handles.imstream(camid).binary_t));

set(handles.sl_area, 'Value', handles.imstream(camid).area_t(1));
set(handles.txt_area_t, 'String', sprintf('Area:%.0f', handles.imstream(camid).area_t(1)));
% 2/2/2011 end change -- retrieve thresholds
[bbox, br]=get_search_box(handles, camid);

if ~isempty(bbox)
    [Z, handles.imstream(camid).seg_img]=getZ(handles.imstream(camid),camid, ...
        handles.optTrack, handles.imstream(camid).binary_t*3/4, bbox, br);
else
    [Z, handles.imstream(camid).seg_img]=getZ(handles.imstream(camid),camid,handles.optTrack);
end

if ~isempty(Z)                                
    % handles.imstream(camid).centroids=cat(1,si(idx).Centroid);
    handles.imstream(camid).centroids=cat(2,Z.u)';
    handles.imstream(camid).Z=Z;
end

% save threshold values
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.bt, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).binary_t;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.at, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).area_t;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.br, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).br;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ns, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).noise_std;




function handles = extract_pts(handles, camid)


[bbox br]=get_search_box(handles, camid);

if ~isempty(bbox)
    [Z handles.imstream(camid).seg_img]=getZ(handles.imstream(camid), ...
        camid,handles.optTrack, handles.imstream(camid).binary_t*3/4, bbox, br);
else
    [Z handles.imstream(camid).seg_img]=getZ(handles.imstream(camid),camid,handles.optTrack);
end

% handles.imstream(camid).centroids=cat(1,si(idx).Centroid);
if ~isempty(Z)                                
    % handles.imstream(camid).centroids=cat(1,si(idx).Centroid);
    handles.imstream(camid).centroids=cat(2,Z.u)';
    handles.imstream(camid).Z=Z;
end
% save threshold values
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.bt, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).binary_t;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.at, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).area_t;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.br, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).br;
% 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ns, 2);
% handles.Z(r,c(camid))=handles.imstream(camid).noise_std;

function [bbox br]=get_search_box(handles, camid)
Xi=handles.Xi;

[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1:Xi.nX, 1);
if handles.backtrack
    if handles.Xh(r(1),c+1)
        zh_=w2cam(handles.Xh(r(Xi.ri),c+1)-handles.Xh(r(Xi.rdi),c+1)*handles.optTrack.dt, handles.get_cam_calib(camid));
        bbox=[max(1,ceil(zh_(1))-handles.searchbbox(1)), ...
              max(1,ceil(zh_(2))-handles.searchbbox(2)), ...
              min(handles.imstream(camid).img_info.Height, handles.searchbbox(1)*2), ...
              min(handles.imstream(camid).img_info.Width,  handles.searchbbox(2)*2)];
        img_plane_vel=norm(handles.Xh(r(Xi.rdi(1:2)),c+1))/1000;
        br=max(1,3-floor(min(3,img_plane_vel)));  
    else
        bbox=[];
        br=[];
    end
else
    if handles.Xh(r(1),c-1)
        zh_=w2cam(handles.Xh(r(Xi.ri),c-1)+handles.Xh(r(Xi.rdi),c-1)*handles.optTrack.dt, handles.get_cam_calib(camid));
        bbox=[max(1,ceil(zh_(1))-handles.searchbbox(1)), ...
              max(1,ceil(zh_(2))-handles.searchbbox(2)), ...
              min(handles.imstream(camid).img_info.Height, handles.searchbbox(1)*2), ...
              min(handles.imstream(camid).img_info.Width,  handles.searchbbox(2)*2)];
        img_plane_vel=norm(handles.Xh(r(Xi.rdi(1:2)),c-1))/1000;
        br=max(1,3-floor(min(3,img_plane_vel)));  
    else
        bbox=[];
        br=[];
    end
end

function show_frames(camid, imgtype, handles)

display_image(camid, imgtype, handles);

% plot markers
plot_curr_mq_k(handles, camid);
if get(handles.cb_see_swarm, 'Value')
    plot_all_mq_k(handles, camid);
end

function display_image(camid, imgtype, handles)

% get the axes
axes(handles.ah(camid));
gca; cla;

if(imgtype==0)
    % read and show the image
    imshow(handles.imstream(camid).orig_img); hold on;
elseif(imgtype==1)
    imshow(handles.imstream(camid).seg_img); hold on;
end

% Get the first click and center around that
set(gca, 'xlim',handles.imstream(camid).xlim);
set(gca, 'ylim',handles.imstream(camid).ylim);


% --- Executes on slider movement.
function sl_area_Callback(hObject, eventdata, handles)
% hObject    handle to sl_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% find max value for each pixel through all images
camid=handles.tunecam;
handles.imstream(camid).area_t(1)=get(hObject,'Value');
set(handles.txt_area_t, 'String', sprintf('Area:%.0f', handles.imstream(camid).area_t(1)));

set(gcf,'Pointer','watch');
handles=extract_pts(handles, camid);

show_frames(camid, get_frame_choice(handles, camid), handles);
% if(camid==1 && get(handles.rb_segmented_L, 'Value'))
%     show_frames(camid, 1, handles);
% elseif(camid==2 && get(handles.rb_segmented_R, 'Value'))
%     show_frames(camid, 1, handles);
% end


% Update handles structure
guidata(hObject, handles);

set(gcf,'Pointer','arrow');

% --- Executes during object creation, after setting all properties.
function sl_area_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_area (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_edge.
function pb_edge_Callback(hObject, eventdata, handles)
% hObject    handle to pb_edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Left Camera
ii=1;
% get the axes
axes(handles.ah(ii));

% read and show the image
Xt=imread([handles.frmloc, handles.imstream(ii).flist(handles.k).name]);

% run canny edge detector
eXt=edge(Xt, 'canny');

Xt=double(Xt)+eXt*10000;
imshow(Xt/65535);

hold on;

% --- Executes on button press in rb_tune_L.
function rb_tune_L_Callback(hObject, eventdata, handles)
% hObject    handle to rb_tune_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
    handles=choose_camera_to_tune(handles,1);
end


guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of rb_tune_L


% --- Executes on button press in rb_tune_R.
function rb_tune_R_Callback(hObject, eventdata, handles)
% hObject    handle to rb_tune_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
    handles=choose_camera_to_tune(handles,2);
end

guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of rb_tune_R


function handles = choose_camera_to_tune(handles, camid)

handles.tunecam=camid;
set(handles.sl_bgt, 'Value', handles.imstream(camid).binary_t);
set(handles.txt_bin_t, 'String', sprintf('Intensity:%.4f', handles.imstream(camid).binary_t));

set(handles.sl_area, 'Value', handles.imstream(camid).area_t(1));
set(handles.txt_area_t, 'String', sprintf('Area:%.0f', handles.imstream(camid).area_t(1)));

% --- Executes on button press in rb_original_R.
function rb_original_R_Callback(hObject, eventdata, handles)
% hObject    handle to rb_original_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
   show_frames(2,0, handles);
end
% Hint: get(hObject,'Value') returns toggle state of rb_original_R


% --- Executes on button press in rb_segmented_R.
function rb_segmented_R_Callback(hObject, eventdata, handles)
% hObject    handle to rb_segmented_R (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
   set(handles.rb_tune_R, 'Value', 1);
   handles=choose_camera_to_tune(handles,2);
   show_frames(2,1, handles);
end
guidata(hObject, handles);


% Hint: get(hObject,'Value') returns toggle state of rb_segmented_R


% --- Executes on button press in rb_segmented_L.
function rb_segmented_L_Callback(hObject, eventdata, handles)
% hObject    handle to rb_segmented_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
   set(handles.rb_tune_L, 'Value', 1);    
   handles=choose_camera_to_tune(handles,1);
   show_frames(1,1, handles);
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of rb_segmented_L


% --- Executes on button press in rb_original_L.
function rb_original_L_Callback(hObject, eventdata, handles)
% hObject    handle to rb_original_L (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(hObject,'Value'))
   show_frames(1,0, handles);
end
% Hint: get(hObject,'Value') returns toggle state of rb_original_L


% --- Executes on button press in pb_cam1_crosshair.
function pb_cam1_crosshair_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cam1_crosshair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

jj=handles.mqid;

axes(handles.ah(1));
% get ready to click on a mosquito
[handles.frame(handles.k).mq(jj).c(1,1), handles.frame(handles.k).mq(jj).c(2,1)]=ginput(1);
plot(handles.frame(handles.k).mq(jj).c(1,1), ...
     handles.frame(handles.k).mq(jj).c(2,1), 'k+','MarkerSize', handles.ms(1));
 
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.uv, 2);
% % handles.Z(r,c(1))=handles.frame(handles.k).mq(jj).c(:,1);
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ua, 2);
% handles.Z(r,c(1))=2;

handles=update(handles);

for c_ii=1:handles.nc
    plot_curr_mq_k(handles, c_ii);
    if get(handles.cb_see_swarm, 'Value')
        plot_all_mq_k(handles, c_ii);
    end
end

guidata(hObject, handles);

% --- Executes on button press in pb_cam2_crosshair.
function pb_cam2_crosshair_Callback(hObject, eventdata, handles)
% hObject    handle to pb_cam2_crosshair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
jj=handles.mqid;
axes(handles.axes2);
handles.frame(handles.k).mq(jj).c(:,2)=ginput(1);
plot(handles.frame(handles.k).mq(jj).c(1,2), ...
     handles.frame(handles.k).mq(jj).c(2,2), 'k+', 'MarkerSize', handles.ms(2));

% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.uv, 2);
% handles.Z(r,c(2))=handles.frame(handles.k).mq(jj).c(:,2);
% [r, c]=getind(handles.Zi.nZ, handles.k, handles.mqid, handles.Zi.ua, 2);
% handles.Z(r,c(2))=2; 

% set_instr('[A] Click on Triangulate to localize the mosquito', handles)
[handles, ind]=t1_triangulate(handles);
if ind
    set_instr('[A] Track spliced to an auto version... Check override or keep moving', handles);
end

[r, c]=getind(handles.Xi.nX, handles.k, jj, 1:handles.Xi.nX, 1);
lsr=handles.Xh(r(handles.Xi.ri), c);

for c_ii=1:handles.nc
    plot_curr_mq_k(handles, c_ii);
    if get(handles.cb_see_swarm, 'Value')
        plot_all_mq_k(handles, c_ii);
    end
end

instr=sprintf('3D position: [%.1f, %.1f, %.1f]\nPress next (>) or save', ceil(lsr(1)), ceil(lsr(2)), ceil(lsr(3)));
set_instr(instr, handles)

guidata(hObject, handles);


% --------------------------------------------------------------------
function mb_pan_Callback(hObject, eventdata, handles)
% hObject    handle to mb_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gca;
pan on;

for ii=1:handles.nc
    axes(handles.ah(ii));
    handles.imstream(ii).xlim=get(gca, 'xlim');
    handles.imstream(ii).ylim=get(gca, 'ylim');
end
guidata(hObject, handles);

% set_instr('Remember to press Normal to remember the magnification in memory', handles);

% --------------------------------------------------------------------
function mb_zoomin_Callback(hObject, eventdata, handles)
% hObject    handle to mb_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gca;
zoom on;

for ii=1:handles.nc
    axes(handles.ah(ii));
    handles.imstream(ii).xlim=get(gca, 'xlim');
    handles.imstream(ii).ylim=get(gca, 'ylim');
end
guidata(hObject, handles);

% set_instr('Remember to press Normal to remember the magnification in memory', handles);

% --------------------------------------------------------------------
function mb_normal_Callback(hObject, eventdata, handles)
% hObject    handle to mb_normal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom off
pan off

% Get axes values for magnification
for ii=1:handles.nc
    axes(handles.ah(ii));
    handles.imstream(ii).xlim=get(gca, 'xlim');
    handles.imstream(ii).ylim=get(gca, 'ylim');
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function uicm_help_Callback(hObject, eventdata, handles)
% hObject    handle to uicm_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% 

% --- Executes on button press in pb_next1.
function pb_next1_Callback(hObject, eventdata, handles)
% hObject    handle to pb_next1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.k+1 > size(handles.Xh,2)-handles.optTrack.br0)
    set_instr('[!] You are on the last possible frame', handles)
else
    handles.k = handles.k + 1;

    set(handles.ed_frame, 'String', sprintf('%d', handles.k));

    handles.backtrack=0;

    % Store all images that are needed to segment
    for ii=1:handles.nc
        handles=store_orig_image(ii, handles);
        % Call function to segment image
        if(get(handles.cb_track, 'Value'))
            handles=update_img_array(handles, ii);
            handles=apply_bg_subtract(handles, ii);
        end
    end


    % Get next best dots from motion model
    % if(handles.imstream(1).stageflags(4,handles.k-1, handles.mqid))
    if(get(handles.cb_track, 'Value'))
        [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1:handles.Xi.nX, 1);
        if(size(handles.Xh,1)>=r(1))
            if(~handles.Xh(r(1),c) || get(handles.cb_override_existing_track, 'Value'))
                handles = predict(handles);
            end
        end
    end
    % end

    % Get axes values for magnification
    for ii=1:handles.nc
        axes(handles.ah(ii));
        handles.imstream(ii).xlim=get(gca, 'xlim');
        handles.imstream(ii).ylim=get(gca, 'ylim');
        display_image(ii,get_frame_choice(handles,ii), handles);
    end
    
    if(get(handles.cb_track, 'Value'))
        if(~handles.Xh(r(1),c)|| get(handles.cb_override_existing_track, 'Value'))
            handles=update(handles);
%         else
%             set_instr('[A] Click on a mosquito in the left frame', handles)
        end
    end

    for c_ii=1:handles.nc
        plot_curr_mq_k(handles, c_ii);
        if get(handles.cb_see_swarm, 'Value')
            plot_all_mq_k(handles, c_ii);
        end
    end
    
    % start change 1/24/2011 : show 3D position and speed
    plot_curr_mq_k_3d_speed(handles);
    % end change 1/24/2011 : show 3D position and speed


    guidata(hObject, handles);

end


function ed_goto_Callback(hObject, eventdata, handles)
% hObject    handle to ed_goto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_goto as text
%        str2double(get(hObject,'String')) returns contents of ed_goto as a double


% --- Executes during object creation, after setting all properties.
function ed_goto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_goto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ed_frame_Callback(hObject, eventdata, handles)
% hObject    handle to ed_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns cont


% --- Executes during object creation, after setting all properties.
function ed_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_track.
function cb_track_Callback(hObject, eventdata, handles)
% hObject    handle to cb_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_track

if(~get(hObject, 'Value'))
    handles=disable_tracking(handles);
else
    handles=enable_tracking(handles);
end

guidata(hObject, handles);    


function handles=disable_tracking(handles)

% handles.k_save=handles.k;
% fprintf('Unchecked tracking at %d...\n', handles.k);
% set(handles.pb_record, 'Enable', 'off');
% set(handles.pb_pair, 'Enable', 'off');
% set(handles.pb_cam1_crosshair, 'Enable', 'off');    
% set(handles.pb_cam2_crosshair, 'Enable', 'off');    
% set(handles.rb_segmented_L, 'Enable', 'off');     
% set(handles.rb_segmented_R, 'Enable', 'off');    
if ~get(handles.rb_original_L, 'Value')
    set(handles.rb_original_L, 'Value', 1);
    show_frames(1, 0, handles);
    plot_curr_mq_k(handles, 1);
end
if ~get(handles.rb_original_R, 'Value')
    set(handles.rb_original_R, 'Value', 1);
    show_frames(2, 0, handles);
    plot_curr_mq_k(handles, 1);
end
    
function handles=enable_tracking(handles)

% handles.k=handles.k_save;
% fprintf('Going back to tracking at %d...\n', handles.k);    
% set(handles.pb_record, 'Enable', 'on');
% set(handles.pb_pair, 'Enable', 'on');
% set(handles.pb_cam1_crosshair, 'Enable', 'on');    
% set(handles.pb_cam2_crosshair, 'Enable', 'on');  
% set(handles.rb_segmented_L, 'Enable', 'on');     
% set(handles.rb_segmented_R, 'Enable', 'on');    
set(handles.ed_frame, 'String', sprintf('%d', handles.k));    

for c_ii=1:handles.nc
    handles=store_orig_image(c_ii, handles);
     handles=init_img_array(handles, c_ii);
    % Call function to segment image
    handles=apply_bg_subtract(handles, c_ii);
    
    show_frames(c_ii,0, handles);
    plot_curr_mq_k(handles, c_ii);
end

% --- Executes on button press in pb_goto_frm.
function pb_goto_frm_Callback(hObject, eventdata, handles)
% hObject    handle to pb_goto_frm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.backtrack=0;
goto_frm = str2double(get(handles.ed_frame, 'String'));

if (goto_frm <=handles.optTrack.br0 || goto_frm >= size(handles.Xh,2)-handles.optTrack.br0)
    set_instr(sprintf('[!] Frame should be between %d - %d', handles.optTrack.br0+1, size(handles.Xh,2)-handles.optTrack.br0), handles);
else
   
    handles.k=goto_frm;

    % find if the current id exists
    [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1, 1);
    
    

    if(handles.Xh(r(1),c))
        fprintf('Found track for current mqid %d...\n', handles.mqid);
    else
        fprintf('Did not find any track for current mqid %d ...\n', handles.mqid);
        [handles.mqid handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);
    end
    
    % Store all images that are needed to segment
    for ii=1:handles.nc
        handles=store_orig_image(ii, handles);
        
        handles.imstream(ii).imgarr=zeros(handles.imstream(ii).img_info.Height, handles.imstream(ii).img_info.Width, ...
                        handles.optTrack.br0*2+1);
        
        handles=init_img_array(handles, ii);
        % Call function to segment image
        handles=apply_bg_subtract(handles, ii);
    end


    [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1:handles.Xi.nX, 1);

    % Get axes values for magnification
    for ii=1:handles.nc
        axes(handles.ah(ii));
        handles.imstream(ii).xlim=get(gca, 'xlim');
        handles.imstream(ii).ylim=get(gca, 'ylim');
    end

    % show next frames
    if (get(handles.rb_original_L, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_L, 'Value'))
        frame_choice=1;
    end
    display_image(1,frame_choice, handles);

    if (get(handles.rb_original_R, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_R, 'Value'))
        frame_choice=1;
    end

    display_image(2,frame_choice, handles);

    for c_ii=1:handles.nc
        plot_curr_mq_k(handles, c_ii);
        if get(handles.cb_see_swarm, 'Value')
            plot_all_mq_k(handles, c_ii);
        end
    end

     % start change 1/24/2011 : show 3D position and speed
    plot_curr_mq_k_3d_speed(handles);
    % end change 1/24/2011 : show 3D position and speed
    
%     set_instr('[I] Click on a mosquito in the left frame', handles)
    guidata(hObject, handles);

end

set(handles.ed_frame, 'String', sprintf('%d', handles.k));

function frame_choice=get_frame_choice(handles, camid)
frame_choice=0;
if camid==1
    if (get(handles.rb_original_L, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_L, 'Value'))
        frame_choice=1;
    end
elseif camid==2
    if (get(handles.rb_original_R, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_R, 'Value'))
        frame_choice=1;
    end
end

% --- Executes on button press in cb_see_swarm.
function cb_see_swarm_Callback(hObject, eventdata, handles)
% hObject    handle to cb_see_swarm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(get(hObject,'Value'))
    for ii=1:handles.nc
        plot_all_mq_k(handles, ii);
    end
else
    for ii=1:handles.nc
        show_frames(ii, 0, handles)
    end
end


% Hint: get(hObject,'Value') returns toggle state of cb_see_swarm


% --- Executes on button press in cb_trajectory.
function cb_trajectory_Callback(hObject, eventdata, handles)
% hObject    handle to cb_trajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


for c_ii=1:handles.nc
    show_frames(c_ii, 0, handles);
    plot_curr_mq_k(handles, c_ii);
end
% Hint: get(hObject,'Value') returns toggle state of cb_trajectory


function keypress_nav(src, event)

% event.Key
% event.Modifier
% event.Character
handles = guidata(src);
switch event.Key
    case 'e'
        pb_next1_Callback(handles.pb_next1, [], handles);    
    case 'q'
        pb_prev1_Callback(handles.pb_prev1, [], handles);    
    case 'space'
        set(handles.pb_forward, 'UserData', 0);
        set(handles.pb_rewind, 'UserData', 0);
    case '1'
        pb_cam1_crosshair_Callback(handles.pb_cam1_crosshair, [], handles);
    case '2'
        pb_cam2_crosshair_Callback(handles.pb_cam2_crosshair, [], handles);
end


% --------------------------------------------------------------------
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cb_override_existing_track.
function cb_override_existing_track_Callback(hObject, eventdata, handles)
% hObject    handle to cb_override_existing_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_override_existing_track


% --- Executes on button press in pb_fplus.
function pb_fplus_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fplus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the first click and center around that
[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, handles.Xi.ri, 1);

if (size(handles.Xh,1) >=r(1) && size(handles.Xh,2)>=c)
    handles.fz=min(8, handles.fz+1);
    for c_ii=1:handles.nc
        handles=focused_zoom(handles, c_ii);
    end

    if (get(handles.rb_original_R, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_R, 'Value'))
        frame_choice=1;
    end

    display_image(2,frame_choice, handles);
    for c_ii=1:handles.nc
        plot_curr_mq_k(handles, c_ii);
        if get(handles.cb_see_swarm, 'Value')
            plot_all_mq_k(handles, c_ii);
        end
    end

    guidata(hObject, handles);
end
    

function handles = focused_zoom(handles, camid)

% Get the first click and center around that
[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, handles.Xi.ri, 1);

if handles.Xh(r,c)
    pix=w2cam(handles.Xh(r,c), handles.get_cam_calib(camid));
    axes(handles.ah(camid));
    box_w=handles.imstream(camid).img_info.Width/handles.fz;
    xmin=max(1,ceil(pix(1)-box_w/2));
    xmax=min(floor(pix(1)+box_w/2), handles.imstream(camid).img_info.Width);
    handles.imstream(camid).xlim=[xmin, xmax];

    box_h=handles.imstream(camid).img_info.Height/handles.fz;
    ymin=max(1,ceil(pix(2)-box_h/2));
    ymax=min(floor(pix(2)+box_h/2), handles.imstream(camid).img_info.Height);
    handles.imstream(camid).ylim=[ymin, ymax];

    set(gca, 'xlim',handles.imstream(camid).xlim);
    set(gca, 'ylim',handles.imstream(camid).ylim);
end


% --------------------------------------------------------------------
function m_help_Callback(hObject, eventdata, handles)
% hObject    handle to m_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function m_browser_Callback(hObject, eventdata, handles)
% hObject    handle to m_browser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

web('https://docs.google.com/leaf?id=0B3kWtXv3eddnYWI4MWM2ZTUtZGE4OS00NWNmLWIwYWUtYjEyOTE2MzJmODEy&hl=en', '-browser')


% --- Executes on button press in pb_fminus.
function pb_fminus_Callback(hObject, eventdata, handles)
% hObject    handle to pb_fminus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the first click and center around that
[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, handles.Xi.ri, 1);

if (size(handles.Xh,1) >=r(1) && size(handles.Xh,2)>=c)
    handles.fz=max(2, handles.fz-1);
    for c_ii=1:handles.nc
        camid=c_ii;

        if handles.Xh(r,c)
            pix=w2cam(handles.Xh(r,c), handles.get_cam_calib(camid));
            axes(handles.ah(camid));
            box_w=handles.imstream(c_ii).img_info.Width/handles.fz;
            xmin=max(1,ceil(pix(1)-box_w/2));
            xmax=min(floor(pix(1)+box_w/2), handles.imstream(c_ii).img_info.Width);
            handles.imstream(camid).xlim=[xmin, xmax];

            box_h=handles.imstream(c_ii).img_info.Height/handles.fz;
            ymin=max(1,ceil(pix(2)-box_h/2));
            ymax=min(floor(pix(2)+box_h/2), handles.imstream(c_ii).img_info.Height);
            handles.imstream(camid).ylim=[ymin, ymax];

            set(gca, 'xlim',handles.imstream(camid).xlim);
            set(gca, 'ylim',handles.imstream(camid).ylim);
        end
    end

    if (get(handles.rb_original_R, 'Value'))
        frame_choice=0;
    elseif(get(handles.rb_segmented_R, 'Value'))
        frame_choice=1;
    end

    display_image(2,frame_choice, handles);
    for c_ii=1:handles.nc
        plot_curr_mq_k(handles, c_ii);
        if get(handles.cb_see_swarm, 'Value')
            plot_all_mq_k(handles, c_ii);
        end
    end

    guidata(hObject, handles);
end


% --- Executes on button press in pb_restore.
function pb_restore_Callback(hObject, eventdata, handles)
% hObject    handle to pb_restore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for c_ii=1:handles.nc
    camid=c_ii;
    
    axes(handles.ah(camid));
    handles.imstream(camid).xlim=[1, handles.imstream(camid).img_info.Width];
    handles.imstream(camid).ylim=[1, handles.imstream(camid).img_info.Height];
    handles.fz=1;
    set(gca, 'xlim',handles.imstream(camid).xlim);
    set(gca, 'ylim',handles.imstream(camid).ylim);
end

guidata(hObject, handles);


% --- Executes on button press in pb_forward.
function pb_forward_Callback(hObject, eventdata, handles)
% hObject    handle to pb_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(handles.pb_forward, 'UserData')
    set(handles.pb_forward, 'UserData', 0);
    set(handles.pb_forward, 'String', '>>');
else
    set(handles.pb_forward, 'UserData', 1)
    set(handles.pb_forward, 'String', 'I I');


    set(handles.cb_track, 'Value', 0);
    handles=disable_tracking(handles);

    % find if the current id exists
    [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1, 1);

    if(handles.Xh(r(1),c))
        % get all the timesteps for the current id
        krange=find(handles.Xh(r(1),:)~=0);
        % call fplus 2-3 times 
        handles.fz=5;
        k=handles.k;
        while k<krange(end) && get(handles.pb_forward, 'UserData')
            handles.k=k;
            for c_ii=1:handles.nc
                handles=store_orig_image(c_ii, handles);
                handles=focused_zoom(handles, c_ii);
                display_image(c_ii, 0, handles)
            end

            for c_ii=1:handles.nc
                plot_curr_mq_k(handles, c_ii);
                if get(handles.cb_see_swarm, 'Value')
                    plot_all_mq_k(handles, c_ii);
                end
            end

            % start change 1/24/2011 : show 3D position and speed
            plot_curr_mq_k_3d_speed(handles);
            k=k+1;
        end
    end

    set(gcf,'Pointer','watch');
    set(handles.cb_track, 'Value', 1);
    set(handles.pb_forward, 'UserData', 0);
    set(handles.pb_forward, 'String', '>>');
    
    handles=enable_tracking(handles);
    guidata(hObject, handles);
end

set(gcf,'Pointer','arrow');

% --- Executes on button press in pb_rewind.
function pb_rewind_Callback(hObject, eventdata, handles)
% hObject    handle to pb_rewind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



if get(handles.pb_rewind, 'UserData')
    set(handles.pb_rewind, 'UserData', 0);
    set(handles.pb_rewind, 'String', '<<');
else
    set(handles.pb_rewind, 'UserData', 1)
    set(handles.pb_rewind, 'String', 'I I');
    set(handles.cb_track, 'Value', 0);
    handles=disable_tracking(handles);

    % find if the current id exists
    [r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1, 1);


    if(handles.Xh(r(1),c))
        % get all the timesteps for the current id
        krange=find(handles.Xh(r(1),:)~=0);
        
        k=handles.k;
            % call fplus 2-3 times 
            handles.fz=5;
        while k>krange(1) && get(handles.pb_rewind, 'UserData')
            handles.k=k;
            for c_ii=1:handles.nc
                handles=store_orig_image(c_ii, handles);
                handles=focused_zoom(handles, c_ii);
                display_image(c_ii, 0, handles)
            end

            for c_ii=1:handles.nc
                plot_curr_mq_k(handles, c_ii);
                if get(handles.cb_see_swarm, 'Value')
                    plot_all_mq_k(handles, c_ii);
                end
            end

            % start change 1/24/2011 : show 3D position and speed
            plot_curr_mq_k_3d_speed(handles);
            k=k-1;
        end

    end
    set(gcf,'Pointer','watch');
    set(handles.pb_rewind, 'UserData', 0);
    set(handles.pb_rewind, 'String', '<<');
    set(handles.cb_track, 'Value', 1);
    handles=enable_tracking(handles);
    guidata(hObject, handles);
    set(gcf,'Pointer','arrow');
end

% --- Executes on button press in pb_chop_tail.
function pb_chop_tail_Callback(hObject, eventdata, handles)
% hObject    handle to pb_chop_tail (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resp=questdlg('Confirm chopping off past points on the track?');

if strcmp(resp, 'Yes')
    if handles.mqid
        handles.Xh=terminate_track(handles.Xh, handles.mqid, handles.k, handles.Xi, 0);
        guidata(hObject, handles);

        for c_ii=1:handles.nc
            show_frames(c_ii, get_frame_choice(handles, c_ii), handles);
        end
    end
end

% --- Executes on button press in pb_jump_forward.
function pb_jump_forward_Callback(hObject, eventdata, handles)
% hObject    handle to pb_jump_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(gcf,'Pointer','watch');
[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1, 1);
if(handles.Xh(r(1),c))
    % get all the timesteps for the current id
    krange=find(handles.Xh(r(1),:)~=0);
    handles.k=krange(end);
    handles=enable_tracking(handles);
    guidata(hObject, handles);
end
set(gcf,'Pointer','arrow');
% --- Executes on button press in pb_jump_back.
function pb_jump_back_Callback(hObject, eventdata, handles)
% hObject    handle to pb_jump_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(gcf,'Pointer','watch');
[r, c]=getind(handles.Xi.nX, handles.k, handles.mqid, 1, 1);
if(handles.Xh(r(1),c))
    % get all the timesteps for the current id
    krange=find(handles.Xh(r(1),:)~=0);
    handles.k=krange(1);
    handles=enable_tracking(handles);
    guidata(hObject, handles);
end
set(gcf,'Pointer','arrow');


% --- Executes on button press in pb_chop_head.
function pb_chop_head_Callback(hObject, eventdata, handles)
% hObject    handle to pb_chop_head (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

resp=questdlg('Confirm chopping off future points on the track?');

if strcmp(resp, 'Yes')
    if handles.mqid
        handles.Xh=terminate_track(handles.Xh, handles.mqid, handles.k, handles.Xi, 1);
        guidata(hObject, handles);

        for c_ii=1:handles.nc
            show_frames(c_ii, get_frame_choice(handles, c_ii), handles);
        end
    end
end

% --------------------------------------------------------------------
function m_new_track_Callback(hObject, eventdata, handles)
% hObject    handle to m_new_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 0);
guidata(hObject, handles);

% --------------------------------------------------------------------
function m_file_Callback(hObject, eventdata, handles)
% hObject    handle to m_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_show_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to m_show_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

show_tracked_mq(handles.Xh, handles.Xi, 0, 1);
% --------------------------------------------------------------------
function m_save_Callback(hObject, eventdata, handles)
% hObject    handle to m_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


logger([handles.floc, '/output/trkrun.log'], sprintf('trackone.m before save k=%d, mqid=%d', handles.k, handles.mqid));

t1_save_tracks(handles)

set_instr('[I] Tracks at current time step saved.', handles);

guidata(hObject, handles);


% --------------------------------------------------------------------
function m_cam2world_Callback(hObject, eventdata, handles)
% hObject    handle to m_cam2world (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fdatfile=[handles.dataloc, sprintf('data_mq_%s_F.mat', handles.expname)];

if exist(fdatfile, 'file')
    cam2world(fdatfile);
else
    fprintf('[I] Did not find smoothed data.. Converting raw data\n');
    cam2world(handles.datfile1);
end

% --------------------------------------------------------------------
function m_filter_and_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to m_filter_and_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


fdatfile=filter_and_smooth(handles.datfile1);


% --------------------------------------------------------------------
function m_exit_Callback(hObject, eventdata, handles)
% hObject    handle to m_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% resp=questdlg('Save before exit?');
% 
% if strcmp(resp, 'Yes')
%     t1_save_tracks(handles.Xh, handles.Z, handles.datfile1, handles.Xi, handles.Zi);
% end
%     
% fprintf('Exiting.....\n');
close(handles.figure1);


% --------------------------------------------------------------------
function m_open_track_k_Callback(hObject, eventdata, handles)
% hObject    handle to m_open_track_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);
% Display image frame in the axis for each camera
for ii=1:handles.nc
    show_frames(ii, 0, handles);
end

plot_curr_mq_k_3d_speed(handles);

guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


resp=questdlg('Save before exit?');

if strcmp(resp, 'Yes')
    t1_save_tracks(handles)
end
% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function ui_new_track_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_new_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 0);
guidata(hObject, handles);


% --------------------------------------------------------------------
function ui_open_track_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_open_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);
% Display image frame in the axis for each camera
for ii=1:handles.nc
    show_frames(ii, 0, handles);
end

plot_curr_mq_k_3d_speed(handles);

guidata(hObject, handles);


% --------------------------------------------------------------------
function ui_save_tracks_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ui_save_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

logger([handles.floc, '/output/trkrun.log'], sprintf('trackone.m before save k=%d, mqid=%d', handles.k, handles.mqid));

t1_save_tracks(handles)

set_instr('[I] Tracks at current time step saved.', handles);

guidata(hObject, handles);


% --------------------------------------------------------------------
function m_delete_track_Callback(hObject, eventdata, handles)
% hObject    handle to m_delete_track (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);

resp=questdlg('Confirm deletion of track?');

if strcmp(resp, 'Yes')
    if handles.mqid
        handles.Xh=terminate_track(handles.Xh, handles.mqid, 1, handles.Xi, 1);
        guidata(hObject, handles);
        for c_ii=1:handles.nc
            show_frames(c_ii, get_frame_choice(handles, c_ii), handles);
        end
    end
end

[handles.mqid, handles.mqids]=t1_get_mqid(handles.Xh, handles.k, handles.Xi, 1);


% --------------------------------------------------------------------
function m_tools_Callback(hObject, eventdata, handles)
% hObject    handle to m_tools (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function m_swap_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to m_swap_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


ids=input('ids of tracks to be swapped [id1 id2], or []=cancel: ');
if ~isempty(ids)
    curr_ids=show_tracked_mq(handles.Xh, handles.Xi, 0, 0);
    if ismember(ids(1), curr_ids) && ismember(ids(2), curr_ids)
        Xi=handles.Xi;
        kswap=input('timestep/frame for swap: ');
        
        [r1, c1]=getind(Xi.nX, kswap, ids(1), 1:Xi.nX, 1);
        [r2, c2]=getind(Xi.nX, kswap, ids(2), 1:Xi.nX, 1);

        tXh=handles.Xh(r2,c2:end);
        handles.Xh(r2,c2:end)=handles.Xh(r1,c1:end);
        handles.Xh(r1,c1:end)=tXh;
        fprintf('[I] ids %d %d swapped after frame %d\n', ids(1), ids(2), kswap);
    else
        fprintf('[!] Ids must belong to a tracked target\n');
    end
end
    

guidata(hObject, handles);
