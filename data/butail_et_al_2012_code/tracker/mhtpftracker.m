function mhtpftracker

close all

run initproc
%reset(RandStream.getDefaultStream)

nt_max=200; 

Xi=strXi;
scan=[];


% Image processing
global IMG
IMG.large=150; % indicating upper threshold
IMG.fg_is_dark=1;
IMG.br0=3;
IMG.min_bt=0.0075;
IMG.epc_t=.5; % epipolar constraint threshold abs(z2'*F*z1)
IMG.swarm_boundaries=swarm_boundaries;

% Filter
global FLT
FLT.lfn=@(x, p) pf_update(x, p, cams, Xi);
FLT.motion=@(x) mq_motion(x, IMG);
FLT.sigma_ep=diag([2 2]); % end points
FLT.Np=10;
FLT.inittarget=@(z) initMq(z, cams, IMG);
FLT.figless=1;
FLT.te=1/40; % exposure time
FLT.dt=1/60;
% half height and width of bounding box for searching
FLT.bbox=[50,50]; 
FLT.fg_is_dark=1;
FLT.br0=5; %

% MHT params
global MHT
MHT.Pd=0.95;
MHT.V=100;
MHT.Bft=1/MHT.V;
MHT.Bnt=1/MHT.V;
MHT.gate=16;
MHT.gatefun=@(x,z,S) (x-z)'/inv(S)*(x-z);
MHT.sortbest=4;
MHT.nscanback=1; % look ns steps back to finalize DA. e.g. ns=3, k=7, k-ns=4
MHT.kbest=0;
MHT.reduction_strategy='bestsort';

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
        imstream(cc).img_info.Width, IMG.br0*2+1);
    bgparams=init_bgparams(imstream(cc), cc, floc, frmloc, ...
                    expected_nmq, IMG);
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

nframes=min(size(imstream(1).flist,1),size(imstream(2).flist,1));

% ------------------------
% Initialize the variables
% ------------------------
for cc=1:nc
    fg{cc}=zeros(imstream(cc).img_info.Height, imstream(cc).img_info.Width);
end

% fg=zeros(imstream(1).img_info.Height, imstream(2).img_info.Width, nc);
k0=IMG.br0+1;
Xh=[];

% the number of particles are number of particles, 
tp=zeros(FLT.Np*(MHT.nscanback+1)*nt_max, Xi.nX); 
tp_=tp;

kF=input(sprintf('Start frame=%d\nEnd frame=? []=%d, #=other: ',k0, nframes-IMG.br0 ));
if isempty(kF)
    kF=nframes-IMG.br0;
end

% ------------------------
% Begin tracker
% ------------------------
for k=k0:kF % for each time step / frame
    fprintf('---------------------------------------------------------------------\n');  
    tic
    
    % ------------------------
    % Upload the images into memory for extracting measurements
    % ------------------------  
    for cc=1:nc % for each camera
        imstream(cc).imgarr=update_imgarr(imstream(cc).imgarr, k, k0, ...
                    imstream(cc).noise_std, imstream(cc).flist, frmloc);
    end
    
    % ---------------------------------------
    % Extract measurements
    % ---------------------------------------
    for cc=1:nc % for each camera
        [imstream(cc).Zk, fg{cc}]=getZ(imstream(cc));
    end

       
    if gdebug, debugger(k, scan, tp_, frmloc, imstream, fg, cams, tc, Xi, FLT.Np, 'pre_update'); end
         
    % --------------------------------------
    % Adaptive thresholding 
    % based on gating we look for mosquitoes that may be thresholded out
    % --------------------------------------
    for cc=1:nc
        imstream(cc)=adaptive_thresholding(tp_, cams(cc), imstream(cc), cc, Xi);
    end
    
    % ---------------------------------------------------
    % Measurement collection is a combination satisfying epipolar
    % constraint (for 2 cameras only)
    % ---------------------------------------------------
    F=get_F_for_stereo(cams);
    nz1=size(imstream(1).Zk,2);
    nz2=size(imstream(2).Zk,2);
    cost=ones(nz1, nz2);
    for z1=1:nz1
        for z2=1:nz2
            cost(z1,z2)=abs([imstream(2).Zk(z2).u; 1]'*F*[imstream(1).Zk(z1).u; 1]);
        end
    end
    idx=find(cost<IMG.epc_t);
    [c1, c2]=ind2sub(size(cost), idx);
    Zk=[imstream(1).Zk(c1)', ...
         imstream(2).Zk(c2)'];
    mfpr(k,sprintf('%d measurement pairs... ', size(Zk,1)));
        
    % ---------------------------------------
    % MHT core algorithm
    % ---------------------------------------
    
    [scan, tp, tp_]=mhtcore(k, k0, scan, tp, tp_, Zk, PF, MHT, Xi);
    
    % ---------------------------------------
    % cleaning up for space
    % ---------------------------------------
    tp=cleanup_space(k, tp, Xi);

    % ---------------------------------------
    % shift
    % ---------------------------------------
    [tp, tp_, Xh]=shiftdown(k, scan, tp, tp_, Xh, MHT, Xi);
    
    
    % ---------------------------------------
    % Predict
    % ---------------------------------------
    tp_=predict1(k, tp, tp_, MHT, PF, Xi);
end

save(datafile, 'Xh', 'Xi', 'cams', 'frmlist');

% dump integrity plot into output/movies
% track_integrity_plot
%     mfpr(k, sprintf('Searching for measurements for target'));
%     ctr=0;
%     icov=zeros(nc,max(ids));
%     for cc=1:nc
%         for t_ii=ids'
%             r=getind(Xi.nX, k, t_ii, 1:Xi.nX, 1);
%             icov(cc, t_ii)=min(geticov(p(r,:), Xi, imstream(cc).Zk, cams(cc), IMG));
%             if icov(cc, t_ii) > IMG.gatesize
%                 fprintf('%d_%d(', t_ii, cc);
%                 imstream(cc)=adaptive_thresholding(p(r, :), cams(cc),...
%                                             imstream(cc), cc, Xi, IMG);
%                 % update the covariance values
%                 icov(cc, t_ii)=min(geticov(p(r,:), Xi, imstream(cc).Zk, cams(cc), IMG));
%                 ctr=ctr+1;
%             end
%         end
%         % empty out the extra measurements
%         imstream(cc).Zk2=[];
%     end

function imstream = adaptive_thresholding(p, cam, imstream, cc, Xi)

global gdebug MHT


bt=imstream.binary_t;
foundit=0;
bbox=[1 1 size(imstream.imgarr,2), size(imstream.imgarr,1)];

% ------- stage 1 create an image with lower bt and list the measurements
if isempty(imstream.Zk2)
    imstream.Zk2=getZ(imstream, cc, 2*bt/3, bbox, imstream.br);
end

% ---------- stage 2 search for the closest measurement and if it's still
% not there then don't do anything (to avoid double counting measurements)
[val, idx]=min(geticov(p, Xi, imstream.Zk2, cam));
% don't change gate size here because it affects in other parts of the code
% extra measurements will be included
if val < MHT.gate
    imstream.Zk=[imstream.Zk, imstream.Zk2(idx)];
    foundit=1;
    if gdebug
        Zplot=cat(2,imstream.Zk2(idx).pixel_list);
        figure(cc+2); gcf; hold on;
        plot(Zplot(:,2), Zplot(:,1), 'c.', 'markersize', 2);
    end
end


if foundit
    fprintf('*)..');
else
    fprintf(' )..');
end


function icov=geticov(p, Xi, Zk, cam)
global FLT

pzh_=w2cam(p(Xi.ri,:), cam);
S=cov(pzh_');
zh_=postest(pzh_,ones(FLT.Np,1),1);
icov=zeros(size(Zk,2),1);
for jj=1:size(Zk,2)
    icov(jj)=(zh_-Zk(jj).u)'/S*(zh_-Zk(jj).u);
end
if det(S) < 10^-5
    %keyboard;
    fprintf('[!] Badly scaled innovation matrix... \n');
    % catching badly scaled matrix
end    


% % alternate testing
% for jj=1:size(Zk,2)
%     wts=p_lfn(Zk(jj), p, cam, Xi, optTrack);
%     pz(jj)=sum(wts)/numel(wts);
% end
% 
% pz=pz/max(pz);
% fprintf('--icov----pz---\n');
% for jj=1:size(Zk,2)
% 	fprintf('%.2f   %.2f\n', icov(jj), pz(jj));
% end

function imgarr=update_imgarr(imgarr, k, k0, noise_std, flist, frmloc)
global IMG
%function imgarr=update_imgarr(imgarr, k, optTrack.br0, noise_std, flist, frmloc)
%
% function updates the imagearr incrementally by reading off images from
% the list of files

% if this is the first image then fill up the imgarr
if (k==k0)
    % fill up all that is possible in the begining and everything in
    % the end
    jj1=1;
    for jj=k-IMG.br0:k+IMG.br0
        img=imread([frmloc, flist(jj).name]);
        if(size(img,3)>1), img=rgb2gray(img); end
        if(noise_std>0), img=filter2(fspecial('gaussian', [3 3], noise_std), img); end
        % index is such that we want k to be right on optTrack.br0+1,
        % therefore we go back optTrack.br0 and add jj and then -1
        imgarr(:,:,jj1)=img;
        jj1=jj1+1;
    end
else
    % do a circshift 
    imgarr=circshift(imgarr, [0,0,-1]);

    % populate the last index with k+optTrack.br0 image
    img=imread([frmloc, flist(k+IMG.br0).name]);
    if(size(img,3)>1), img=rgb2gray(img); end
    if(noise_std>0), img=filter2(fspecial('gaussian', [3 3], noise_std), img); end        
    imgarr(:,:,2*optTrack.br0+1)=img;
end


function [Z, fg]=getZ(imstream, varargin)

global gdebug
global IMG

[h, w]=size(imstream.imgarr(:,:,1));

if nargin>1
    imstream.binary_t=varargin{2};
    
    % --------- we can select a specific region instead
    bbox=varargin{3};
    umin=max(bbox(2), imstream.roi(2)); % rows
    umax=min(bbox(2)+bbox(4), imstream.roi(4));
    vmin=max(bbox(1), imstream.roi(1)); % columns
    vmax=min(bbox(1)+bbox(3), imstream.roi(3));
    
    % for a slow mosquito br is large
    % for a fast mosquito br is small
    imstream.br=varargin{4};
    
    imstream.imgarr=imstream.imgarr(umin:umax,vmin:vmax,:);
end

if(IMG.fg_is_dark)
    bg=max(imstream.imgarr(:,:,IMG.br0+1-imstream.br:IMG.br0+1+imstream.br), [], 3);
else
    bg=min(imstream.imgarr(:,:,IMG.br0+1-imstream.br:IMG.br0+1+imstream.br), [], 3);
end

fg=imsubtract(bg, imstream.imgarr(:,:,IMG.br0+1));
fg=fg/imstream.bitval;

fg=im2bw(fg, imstream.binary_t);
if strcmp(version('-release'), '2007a')
    Li=bwlabel(fg);
else
    Li=logical(fg);
end

if nargin ==3
    fg=setRoi(fg, imstream.roi);
end
si = regionprops(Li,'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                    'Orientation', 'PixelIdxList');


%--------------- speeding it up (calling regionprops only once)...
area = [si.Area];
idx = find(area >= imstream.area_t(1) & area < imstream.area_t(2));
new_fg=zeros*fg;
new_fg(cat(1,si(idx).PixelIdxList))=1;

centroids=cat(1,si(idx).Centroid);

if nargin>3 && ~isempty(centroids)
    centroids(:,1)=centroids(:,1)+vmin;
    centroids(:,2)=centroids(:,2)+umin;
end

majora=cat(1,si(idx).MajorAxisLength);
minora=cat(1,si(idx).MinorAxisLength);
orientation=cat(1,si(idx).Orientation);
area=cat(1,si(idx).Area);

[hi, wi]=size(imstream.imgarr(:,:,1));

if nargin>3
    tmp_fg=zeros(h,w);
    tmp_fg(umin:umax,vmin:vmax)=new_fg;
    fg=tmp_fg;
else
    fg=new_fg;
end


if  gdebug % for debug 
    figure(10); gcf;
    j=(uint16(imstream.imgarr(:,:,IMG.br0+1))); %histeq
%     imshow(j); hold on;
    imshow((1-logical(fg))*255); hold on;
    
end

Z=strZ(numel(idx));
for ii=1:numel(idx)
    Z(ii).u=centroids(ii,:)';
    % ------- orientation vector
%     jj=[majora(ii)*cosd(orientation(ii)); -majora(ii)*sind(orientation(ii))]/2;
%     Z(ii).ep=[Z(ii).u+jj Z(ii).u-jj];
    Z(ii).v=[majora(ii).*cosd(orientation(ii)), -majora(ii).*sind(orientation(ii))]';
    % ------- frame transformation bTw = [cost +sint; -sint cost],
    % therefore wTb = bTw'
    Z(ii).sigma=sqrt(abs([majora(ii).*cosd(orientation(ii)) - minora(ii).*sind(orientation(ii)); 
                          majora(ii).*sind(orientation(ii)) + minora(ii).*cosd(orientation(ii))]));
    Z(ii).area=area(ii);
    [i1, i2]=ind2sub([hi, wi], si(idx(ii)).PixelIdxList);

    if nargin>3
        Z(ii).pixel_list=[i1+umin-1, i2+vmin-1];
    else
        Z(ii).pixel_list=[i1,i2];
    end
    
    Z(ii)=setEndPointVelocities(Z(ii));
    % using the index of position as the id
    Z(ii).id=sub2ind([h,w], ceil(Z(ii).u(2)), ceil(Z(ii).u(1)));
end


function [scan, tP, tP_]=mhtcore(k, k0, scan, tP, tP_, Zk, PF, MHT, Xi)

global gf;
global idl;

% ---------------------------------------
% Generate hypotheses
% ---------------------------------------
% try
scan(k).hyp=gen_new_hyp(k, k0, scan, tP_, Zk, PF, MHT, Xi);
% catch
%     k
%     keyboard
% end
% ---------------------------------------
% Compute the probabilities 
% ---------------------------------------

if k>k0 && ~isempty(scan(k-1).hyp)
    for hh=1:size(scan(k).hyp,2)
        pid=scan(k).hyp(hh).parent;

        scan(k).hyp(hh)=hypothesis_prob(k, scan(k).hyp(hh), ...
                            scan(k-1).hyp(pid), tP_, Zk, MHT, PF, Xi);
    end
end
scan = normalize_prob(k, scan);

% ---------------------------------------
% Hypotheses reduction: sort, scanback
% ---------------------------------------
scan = hypotheses_reduction(k, scan, MHT);
scan = nscanback(k, scan, MHT); % scanback

hypmat=cat(1,scan(k).hyp.assignment);
hypmat=assign_new_ids(hypmat);

for hh=1:size(scan(k).hyp,2)
    scan(k).hyp(hh).assignment=hypmat(hh,:); 
end


% ---------------------------------------
% Update hypotheses from previous timestep according to current
% assignments
% ---------------------------------------
tke=[];
tk=unique(hypmat); tk=tk(tk~=0);

if k > k0 && ~isempty(scan(k-1).hyp)
    try
    akp=cat(1,scan(k-1).hyp.assignment);
    catch
        keyboard
    end
    tke=tk(ismember(tk,unique(akp)));
    if size(tke,1)>1, tke=tke'; end
end
if size(tk,1)>1, tk=tk'; end
tkn=tk(~ismember(tk,tke));

for t=tke
    if ~isempty(t)
    [~, zt]=ind2sub(size(hypmat), find(hypmat==t));
    zt=unique(zt);
    [tXh, tP, hypmat]=update1(k, hypmat, zt, t, Zk, tXh, tP, tXh_, tP_, MHT, PF, Xi);
    end
end

for t=tkn
    [~, zt]=ind2sub(size(hypmat), find(hypmat==t));
    zt=unique(zt);
    try
    [tXh, tP]=initialize1(k, zt, t, Zk, tXh, tP, tXh_, tP_, MHT, PF, Xi);
    catch
        keyboard
    end
end

% update hypmat
for hh=1:size(scan(k).hyp,2)
    scan(k).hyp(hh).assignment=hypmat(hh,:);
end

function hyp = gen_new_hyp(k, k0, scan, tp_, Zk, Filter, MHT, Xi)
global gf;
global idl;

k1=min(k, MHT.nscanback+1); % to search the temporary state space

if ~isempty(scan) % existing targets
    hypkm1=scan(k-1).hyp;
    hyp=[];
    for hm1=1:size(hypkm1,2)
        assignment=hypkm1(hm1).assignment;

        tids=assignment(1,assignment(1,:)~=0);

        nz=size(Zk,1);


        % gating
        gate=ones(nz, numel(tids))*MHT.gate*2;
        S=zeros(2,2,numel(tids));
        for t=1:numel(tids)
            [tr, br]=id2tbr(tids(t),gf);
            
            [r, c]=getind(Xi.nX, k1, tr, 1:Xi.nX, Xi.nX);
            
            S(:,:,t)=Filter.H*P_(r,c,br)*Filter.H'+Filter.R;

            for z=1:nz
                gate(z,t)=MHT.gatefun((Xh_(r(Xi.ri),k1,br)), Zk(z,:)', S(:,:,t));
            end
        end

        if MHT.kbest || size(gate,1) > 5
            fprintf('\n #z=%d, running murtys kbest', size(gate,1));
            [hypmat, cost]=murtykbest(gate, @munkres, 4);
            for jj=1:size(gate,2)
                hypmat(hypmat==jj)=tids(jj);
            end
            hyp=[hyp strhyp(hypmat, hm1)];
        else
            % validation matrix
            vm=[ones(nz,1) genvmat(gate, MHT.gate)];
            % hypothesis
            currhyp=strhyp(hypgen2(vm, tids), hm1);
            hyp=[hyp currhyp];
        end

    end
else
   nz=size(Zk,1);
   hyp=get_newid(nz);
end

function hc=hypothesis_prob(k, hc, hp, Xh_, P_, Zk, MHT, Kalman, Xi)

global gf;
if ~hc.prob
hca=hc.assignment;
hpa=hp.assignment;
   
k1=min(k, MHT.nscanback+1); % to search the temporary state space

%target indices of those targets that are existing
t_DT=hca(ismember(hca, hpa));
t_DT=t_DT(t_DT~=0);

% Number of measurements associated with existing targets 
N_DT=numel(t_DT);


% Number of measurements associated false targets
N_FT=sum(hca==0);

% Number of measurements associated with new targets
N_NT=numel(hca)-N_DT-N_FT;

% number of existing targets from previous hyp
N_TGT=sum(hpa>0);

% Equation 16 in paper computing probability of existing targets
prob_DT=1;
for t=1:N_DT
    tid=t_DT(t);
    [tr, br]=id2tbr(tid,gf);
    [r, c]=getind(Xi.nX, k1, tr, 1:Xi.nX, Xi.nX);
    S=Kalman.H*P_(r,c,br)*Kalman.H'+Kalman.R;
    % hca==tid is the measurement id 
    prob_DT=prob_DT*normal2(Zk(hca==tid,:)'-Kalman.H*Xh_(r,k1,br), 0, S);
end

% equation 16 in paper
hc.prob=MHT.Pd^N_DT*(1-MHT.Pd)^(N_TGT-N_DT)*MHT.Bft^N_FT*...
                MHT.Bnt^N_NT*prob_DT*hp.prob;

end


function [p_, wts]=pf_update(Z, p_, cams, Xi)

Np=size(p_,1);
wts=ones(1, Np);
for cc=1:size(cams,2)
    wts=wts.*p_lfn(Z(cc), p_, cams(cc), Xi);
end
% RBPF style approx
Uk=cat(2, Z.u);
if Uk(1,1) && Uk(1,2)
    r=lsTriangulate(Uk, cams);
    p_(Xi.ri,:)=r*ones(1,Np);
    
    wts=wts/sum(wts);
    neff=1/sum(wts.^2);
    if (neff <=  Np/2)
        p_=p_(:,resample(wts));
    end
end

function wts=p_lfn(Zk, p_, cam, Xi)

% the standard deviation in end points is based on location of an end-point
% pixel of a fixed grain of rice on a pendulum in a plane parallel to the camera plane
% that is filmed at 25 frames per second. Since we know the length of the
% grain of rice to be fixed, we assume that any variation in length is due
% to noise in end-point computation. Therefore assuming Gaussian noise,
% noise(length) = noise(e1) + noise(e2). Assuming both are the same,
% std(length) = sqrt(std(ep))
%
% The noise is center position is a function of the streak length. Since
% the bounding ellipse also represents a normal distribution around the
% center of the streak, we can equate the major and minor axis length as
% diag([1/a^2, 1/b^2]) = Sigma^-1. Which implies sigma_x=a, sigma_y=b; BUT
% since the ellipse is not necessarily in the image plane frame we perform
% a transformation [cos(t) -si

wts=p_pos(Zk.u, p_(Xi.ri,:), cam, diag(Zk.sigma)).* ...
    p_mq_velocity(Zk, p_, cam, Xi).*...
    pdf('unif', sqrt(sum(p_(Xi.rdi,:).^2)), 100, 4000);

%     
%               diag(abs(imstream(cc).Zk(assoc_t(t_ii,cc)).v)/2)).*...
%     pdf('unif', p_(r(Xi.ri(3)),:), optTrack.auto.swarm_boundaries(1), optTrack.auto.swarm_boundaries(2));

function wts = p_mq_velocity(Z, tp, cam, Xi)
%function wts = p_mq_velocity(Z, tp, cam, Xi, optTrack)
global FLT

sigma_ep=FLT.sigma_ep;
te=FLT.te;

r=tp(Xi.ri,:);
rdot=tp(Xi.rdi,:);%.*(ones(3,1)*tp(Xi.s,:));
% rddot=tp(Xi.rddi,:); 

% r1=r-rdot*te/2-rddot*te^2/4;
% r2=r+rdot*te/2+rddot*te^2/4;

r1=r-rdot*te/2;
r2=r+rdot*te/2;

e1=w2cam(r1,cam);
e2=w2cam(r2,cam);


wts=normal2(e1,Z.ep(:,1), sigma_ep)'.*normal2(e2,Z.ep(:,2), sigma_ep)' + ...
    normal2(e1,Z.ep(:,2), sigma_ep)'.*normal2(e2,Z.ep(:,1), sigma_ep)';


% wts=normpdf(e1(1,:), Z.h(1), sigma_ep(1)).* normpdf(e1(2,:), Z.h(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.t(1), sigma_ep(1)).* normpdf(e2(2,:), Z.t(2), sigma_ep(2));

% wts=normpdf(e1(1,:), Z.h(1), sigma_ep(1)).* normpdf(e1(2,:), Z.h(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.t(1), sigma_ep(1)).* normpdf(e2(2,:), Z.t(2), sigma_ep(2)) + ...
%     normpdf(e1(1,:), Z.t(1), sigma_ep(1)).* normpdf(e1(2,:), Z.t(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.h(1), sigma_ep(1)).* normpdf(e2(2,:), Z.h(2), sigma_ep(2));
