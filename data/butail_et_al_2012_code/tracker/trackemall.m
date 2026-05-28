function trackemall(floc, frmloc)
% function trackemall(floc, frmloc)
%
% floc is the location of the experiment directory
% frmloc is where the frames are located, this may be different from
% floc/frames
%
% e.g.
% trackemall('../Dropbox/20160326_190755/', '~/data/20160326_190755/frames/');

if nargin < 1
    [floc, frmloc]=getexpdirs;
end

[optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams]=initproc(floc, frmloc);

mainfcn(optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams, floc, frmloc);

% main functin
function mainfcn(optTrack, swarm_boundaries, camid_suffixes, ...
 Xi, dataloc, movieloc, calibloc, expname, ...
 image_id, nc, imageIds, get_cam_calib, cams, floc, frmloc)

global gdebug

% reset(RandStream.getDefaultStream)


% Autotrack options
optTrack.auto.Np=200;
optTrack.auto.large=150; % indicating upper threshold
optTrack.auto.pfout=2; % choices are 1=mean, 2=map
optTrack.auto.max_tt=500;
optTrack.auto.max_ct=500;
optTrack.auto.tt_tl=3; % tentative target track length
optTrack.auto.gatesize=16; % pp 152 Bar-Shalom
optTrack.auto.min_bt=0.0075;
optTrack.auto.epc_t=1; % epipolar constraint threshold abs(z2'*F*z1)
optTrack.auto.da='lgnn'; % options are nnda, gnn, jpda, mht
optTrack.auto.ri=[5 5 10]';
optTrack.auto.rdi=10;
optTrack.auto.anstd=1;
optTrack.auto.mm='cv'; 
optTrack.auto.dt=1/60;
optTrack.auto.sigma_ep=diag([2 2]); % end points
optTrack.auto.swarm_boundaries=swarm_boundaries;


datafile=[dataloc, sprintf('data_mq_auto_%s.mat', expname)];

% initialize the file lists ** store these so that they are checked every
% time any part of the tracker is run
if ~exist(sprintf('%scalib/cam1_bgparams.csv', floc), 'file')
    expected_nmq=input('Max # of mosquitoes you expect to see in this swarm (approx.): ');
else
    expected_nmq=40;
end

imstream(1:nc)=struct('flist', [], 'img_info', [], 'bitval', 0, 'imgarr', [], 'Zk', [], ...
                      'pix_list', [], 'unassoc_Zk', [], 'binary_t', 0, 'area_t', 0, ...
                      'br',0, 'noise_std', 0, 'Zk2', []);
for cc=1:nc
    cam_id=[image_id, camid_suffixes(cc)];
    imstream(cc).flist=dir([frmloc, cam_id,'*.*']);
    imstream(cc).img_info=imfinfo([frmloc, imstream(cc).flist(1).name]);
    %******MATLAB doesn't read the bitdepth properly ??%
    if imstream(cc).img_info.BitDepth >16, imstream(cc).img_info.BitDepth=8; end
    imstream(cc).bitval=2^imstream(cc).img_info.BitDepth-1;
    imstream(cc).imgarr=zeros(imstream(cc).img_info.Height, ...
        imstream(cc).img_info.Width, optTrack.br0*2+1);
    bgparams=init_bgparams(imstream(cc), cc, floc, frmloc, ...
                    expected_nmq, optTrack);
    imstream(cc).binary_t=bgparams.binary_t;
    imstream(cc).area_t=bgparams.area_t;
    imstream(cc).br=bgparams.br;
    imstream(cc).noise_std=bgparams.noise_std;
    imstream(cc).roi=bgparams.roi;
end

try
	frmlist=cat(2,imstream.flist);
catch ME
	fprintf('[!] number of frames from both cameras should be same\n');
	return;
end
% colormap
tc=rand(optTrack.auto.max_ct,3); close all; 
nframes=min(size(imstream(1).flist,1),size(imstream(2).flist,1));

% ------------------------
% Initialize the variables
% ------------------------
for cc=1:nc
    fg{cc}=zeros(imstream(cc).img_info.Height, imstream(cc).img_info.Width);
end
% fg=zeros(imstream(1).img_info.Height, imstream(2).img_info.Width, nc);
Xh=zeros(Xi.cX(1)*optTrack.auto.max_ct, nframes);
p=zeros(Xi.cX(1)*optTrack.auto.max_ct, optTrack.auto.Np);
k0=optTrack.br0+1;

kF=input(sprintf('Start frame=%d\nEnd frame=? []=%d, #=other: ',k0, nframes-optTrack.br0 ));
if isempty(kF)
    kF=nframes-optTrack.br0;
end

logger([floc, '/output/trkrun.log'], sprintf('trackemall.m, start, k0=%d, kF=%d', k0, kF));
logger([floc, '/output/trkrun.log'], sprintf('trackemall.m, start, floc=%s, image_id=%s', floc, image_id));


% ------------------------
% Begin tracker
% ------------------------
try
    for k=k0:kF % for each time step / frame
        fprintf('---------------------------------------------------------------------\n');  

        if k>9999, gdebug=1; end

        % ------------------------
        % Upload the images into memory for extracting measurements
        % ------------------------
        for cc=1:nc % for each camera
            imstream(cc).imgarr=update_imgarr(imstream(cc).imgarr, k, k0, ...
                        imstream(cc).noise_std, imstream(cc).flist, frmloc, optTrack);
        end

        % ------------------------------------------------
        % Extract measurements at current time step
        % ------------------------------------------------ 
        for cc=1:nc % for each camera
            [imstream(cc).Zk, fg{cc}]=getZ(imstream(cc), cc, optTrack);
%             [imstream(cc).Zk, fg(:,:,cc)]=getZ(imstream(cc), cc, optTrack);
        end

        % ------------------------  Get ids at this time 
        ids=getids(Xh, k, Xi);
        mfpr(k, sprintf('total targets=%d', numel(ids)));


        if gdebug, debugger(Xh, Xi, k, p, frmloc, imstream, fg, cams, tc, 'pre_update'); end


        % ------------------------------------------------
        % Adaptive thresholding 
        % based on gating we look for mosquitoes that may be thresholded out
        % ------------------------------------------------
        mfpr(k, sprintf('Searching for measurements for target'));
        ctr=0;
        icov=zeros(nc,max(ids));
        for cc=1:nc
            for t_ii=ids'
                r=getind(Xi.nX, k, t_ii, 1:Xi.nX, 1);
                icov(cc, t_ii)=min(geticov(p(r,:), Xi, imstream(cc).Zk, cams(cc), optTrack));
                if icov(cc, t_ii) > optTrack.auto.gatesize
                    fprintf('%d_%d(', t_ii, cc);
                    imstream(cc)=adaptive_thresholding(p(r, :), cams(cc),...
                                                imstream(cc), cc, Xi, optTrack);
                    % update the covariance values
                    icov(cc, t_ii)=min(geticov(p(r,:), Xi, imstream(cc).Zk, cams(cc), optTrack));
                    ctr=ctr+1;
                end
            end
            % empty out the extra measurements
            imstream(cc).Zk2=[];
        end

        if gdebug
            fprintf('\ninnovation covariance target(cam1, cam2)....\n');
            for t_ii=ids'
                fprintf('%.2d(%.2f, %.2f)\n', t_ii, icov(1,t_ii), icov(2,t_ii));
            end
        end
        if ctr, fprintf('\n'); end


        % ------------------------------------------------
        % Measurement collection is a combination satisfying epipolar
        % constraint (for 2 cameras only)
        % ------------------------------------------------
        F=get_F_for_stereo(cams);
        nz1=size(imstream(1).Zk,2);
        nz2=size(imstream(2).Zk,2);
        cost=ones(nz1, nz2);
        for z1=1:nz1
            for z2=1:nz2
                cost(z1,z2)=abs([imstream(2).Zk(z2).u; 1]'*F*[imstream(1).Zk(z1).u; 1]);
            end
        end
        idx=find(cost<optTrack.auto.epc_t);
        [c1, c2]=ind2sub(size(cost), idx);
        Z3d=[imstream(1).Zk(c1)', ...
             imstream(2).Zk(c2)'];
        mfpr(k,sprintf('%d measurement pairs... ', size(Z3d,1)));


        % ------------------------------------------------ 
        % Data association + Occlusion + updates
        % NOTE: occlusions are detected based on same measurement being
        % assigned to two targets. 
        % ------------------------------------------------
        [allwts, Xh, p, unassigned]=mqda_custom(Xh,p, Z3d, cams ,k, optTrack);


        if gdebug, debugger(Xh, Xi, k, p, frmloc, imstream, fg, cams, tc,  'post_update'); end

        % ------------------------------------------------
        % Remove duplicate tracks ** this is heuristic and should be done away
        % with
        % ------------------------------------------------
    %     [Xh, p, ntem]=terminate_duplicate_tracks(Xh, p, k, optTrack);
    %     mfpr(k,sprintf('%d targets terminated [duplicate]... ', ntem));


        % ------------------------ Reget the ids
        ids=getids(Xh,k,Xi);

        % ------------------------
        % Motion model
        % ------------------------
        for t_ii=ids'
            r=getind(Xi.cX(1), k, t_ii, 1:Xi.cX(1), 1);
            if p(r(Xi.ri(3)),1)
                p(r,:)=mq_motion(p(r,:), optTrack.auto.mm, optTrack.auto.dt);    
                % Using the particle with highest weight
                Xh(r,k+1)=postest(p(r,:), allwts(t_ii,:), optTrack.auto.pfout);        
            end
        end


        % ------------------------
        % New targets
        % ------------------------
        [Xh, p, new_t]= findNewTargets(Xh, p, k+1, cams, Z3d, unassigned, optTrack);
        mfpr(k, sprintf('%d new targets found ...', new_t));
        if gdebug, debugger(Xh, Xi, k, p, frmloc, imstream, fg, cams, tc, 'new_targets'); end
    end
catch ME
    fprintf('[!] check error.log in %s\n', floc);
    fid=fopen([floc, '/output/error.log'], 'a'); 
    fprintf(fid, '\n[!] script:%s\n error:%s, line %d\n', ...
                ME.stack(1).file, ME.message, ME.stack(1).line); 
    fclose(fid); 
    save_data(datafile, Xh, p, Xi, optTrack, cams, frmlist);
end

logger([floc, '/output/trkrun.log'], sprintf('trackemall.m, end saving ...'));
save_data(datafile, Xh, p, Xi, optTrack, cams, frmlist);

function save_data(datafile, Xh, p, Xi, optTrack, cams, frmlist)

save(datafile, 'Xh', 'p', 'Xi', 'optTrack', 'cams', 'frmlist');

% dump integrity plot into output/movies
% track_integrity_plot

