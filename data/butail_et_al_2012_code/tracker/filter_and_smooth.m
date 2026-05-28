function fdatfile=filter_and_smooth(varargin)
% This function uses a constant velocity model on top of the 3D estimates to smooth the velocity
% estimates.

optTrack=config;

% design params
dt=optTrack.dt;
vd=300;
v0=500;

% error in r_1,r_2,r_3 directions (r_3 is along camera optical axis so
% tends to be most erroneous)
Rk=diag([2 2 4]);  

if nargin==1
    % dataloc
%     pathname=varargin{1};
    % data_mq_XXX.mat
    filename=varargin{1};
    [pathname, name]=fileparts(filename);
    
    load(filename, 'Xh', 'Xi');
end

if (strcmp(name(end-1:end), '_F') || strcmp(name(end-1:end), '_W'))
    fprintf('[!] Smooth raw data only...\n');
    return;
end


Xh_Z=Xh; % because Z is Xh only

a=textscan(pathname, '%s', 'Delimiter', '/');
bb=a{:};
floc=sprintf('%s/', bb{1:end-2});


% Processing
% Calibration information
addpath([floc '/calib']);
calibfun=dir([floc 'calib/' '*calib*.m']);
if(size(calibfun,1))
    get_cam_calib=str2func(calibfun.name(1:end-2));
else
    fprintf('[?] Could not find calibration for the cameras. Is this 2D tracking? ...\n');
end


%% 
m_flags=Xh(Xi.ri(3):Xi.nX:end,:);
m_ind=find(m_flags~=0);
[m_mqid, k]=ind2sub(size(m_flags), m_ind);
m_mqid=unique(m_mqid);

fprintf('Smoothing the estimates...\n');

% kalman filter matrices
Fk=  [ 1, 0, 0, dt,  0,  0
              0, 1, 0,  0, dt,  0
              0, 0, 1,  0,  0, dt
              0, 0, 0,  1,  0,  0
              0, 0, 0,  0,  1,  0
              0, 0, 0,  0,  0,  1];

Qk=  [ (dt^3*vd^2)/3,             0,             0, (dt^2*vd^2)/2,             0,             0
                          0, (dt^3*vd^2)/3,             0,             0, (dt^2*vd^2)/2,             0
                          0,             0, (dt^3*vd^2)/3,             0,             0, (dt^2*vd^2)/2
              (dt^2*vd^2)/2,             0,             0,       dt*vd^2,             0,             0
                          0, (dt^2*vd^2)/2,             0,             0,       dt*vd^2,             0
                          0,             0, (dt^2*vd^2)/2,             0,             0,       dt*vd^2]; 
                     

nr1=3;
            
Hk=[eye(nr1), zeros(nr1)];


Klmn.F = Fk;
Klmn.Q= Qk;
Klmn.H = Hk;
Klmn.R= Rk;

RTS=Klmn;



fprintf('Listing mosquito ids found\n');
for jj=m_mqid'
    fprintf('%d ...', jj);
    if ~mod(jj,10)
        fprintf('\n');
    end
end

lcam=get_cam_calib(1);
rcam=get_cam_calib(2);

for mq_ii=m_mqid'
%     fprintf('\n[%.3d]\tr_1\tr_2\tr_3\t|Lcam|\t|Rcam|\t\n===================\n', mq_ii);
    [r, c]=getind(Xi.nX, 1, mq_ii, 1:Xi.nX, 1);
   
    nz_k=find(Xh(r(1),:)~=0);

    k0=nz_k(1);
    kF=nz_k(end);
    
    if kF-k0+1 ~= numel(nz_k)
        fprintf('[!] Incomplete sequence for mqid=%d...\n', mq_ii);
        chk=diff(nz_k,1,2);
        missing_k=find(abs(chk)~=1);
        
        fprintf('Check [%d] frame (no 3D estimate found)\n',nz_k(missing_k));
        cont=input('Continue (and average) []=yes, other=no?: ');
        if ~isempty(cont)
            return
        else
            for ii=1:6
                Xh(r(ii),:)=sma(Xh(r(ii),:), 2);
                Xh(r(ii),isnan(Xh(r(ii),:)))=0;
            end
        end
    end
    
    % covariance
    P_(:,:,k0, mq_ii)=eye(2*nr1)*.5;

    % kalman filter
    for k=k0:kF
        
        % initialize
        [r, c]=getind(Xi.nX, k, mq_ii, 1:Xi.nX, 1);        
        if (k ==k0)
            % r
            Xh_(r(Xi.ri),c)=Xh(r(Xi.ri),c);
            
            [r_p1, c_p1]=getind(Xi.nX, k+1, mq_ii, 1:Xi.nX, 1);
            % r dot
            if Xh(r_p1(1))
                Xh_(r(Xi.rdi),c)=(Xh(r_p1(Xi.ri),c_p1)-Xh(r(Xi.ri),c));
            else
                Xh_(r(Xi.rdi),c)= ones(3,1)*v0;
            end
        end
        
        % update
        Z=Xh_Z(r(Xi.ri),c);
        [Xh(r([Xi.ri, Xi.rdi]),c), P(:,:,k,mq_ii)]=kalmanUpdate(Xh_(r([Xi.ri, Xi.rdi]),c), ...
                                                        P_(:,:,k,mq_ii), Z, Klmn);
                                                    
        lcam_pix_diff(mq_ii).v(k)=norm(w2cam(Xh(r(Xi.ri),c),lcam)-w2cam(Z,lcam));
        rcam_pix_diff(mq_ii).v(k)=norm(w2cam(Xh(r(Xi.ri),c),rcam)-w2cam(Z,rcam));
                                                    
        % predict
        [r1, c1]=getind(Xi.nX, k+1, mq_ii, 1:Xi.nX, 1);
        [Xh_(r1([Xi.ri, Xi.rdi]),c1), P_(:,:,k+1,mq_ii)]= kalmanPredict(Xh(r([Xi.ri, Xi.rdi]),c), ...
                                                                P(:,:,k,mq_ii), Klmn);
    end
    
  
    % Rauch, H. E., Tung, F., & Striebel, C. T. (1965). 
    % Maximum likelihood estimates of linear dynamic systems. 
    % AIAA journal, 3(8), 1445-1450.
    % N here means the last step
   
    for k=kF-1:-1:k0
        
        [r, c]=getind(Xi.nX, k, mq_ii, 1:Xi.nX, 1); 
        Z=Xh_Z(r(Xi.ri),c);
        
        Ck=P(:,:,k,mq_ii)*RTS.F'/P_(:,:,k+1,mq_ii);
        Xh(r([Xi.ri, Xi.rdi]),c)=Xh(r([Xi.ri, Xi.rdi]),c) + ...
            Ck*(Xh(r([Xi.ri, Xi.rdi]),c+1) - RTS.F*Xh(r([Xi.ri, Xi.rdi]),c));
        P(:,:,k,mq_ii)=P(:,:,k,mq_ii) + ...
            Ck*(P(:,:,k+1,mq_ii) - P_(:,:,k+1,mq_ii))*Ck';
        
        lcam_pix_diff(mq_ii).v(k)=norm(w2cam(Xh(r(Xi.ri),c),lcam)-w2cam(Z,lcam));
        rcam_pix_diff(mq_ii).v(k)=norm(w2cam(Xh(r(Xi.ri),c),rcam)-w2cam(Z,rcam));
        
    end

    
        
    
end

figure(1); gcf; clf;
mqc=colormap(lines(max(m_mqid)));

for mq_ii=m_mqid'
    subplot(2,1,1); gca;
    vals=lcam_pix_diff(mq_ii).v;vals=vals(vals~=0);
    plot(vals, 'Color', mqc(mq_ii,:));
    hold on;  
   
    
    subplot(2,1,2); gca;
    vals=rcam_pix_diff(mq_ii).v;vals=vals(vals~=0);
    plot(vals, 'Color', mqc(mq_ii,:));
    hold on;
    
end

subplot(2,1,1); gca; set(gca, 'ylim', [0, 50]);
grid on;
set(gca, 'fontsize', 16);
ylabel('Pixel error (Left camera)');
legend({num2str(m_mqid)});

subplot(2,1,2); gca; set(gca, 'ylim', [0, 50]);
grid on;
set(gca, 'fontsize', 16);
xlabel('timestep(frame)');
ylabel('Pixel error (Right camera)');


saveas(1, sprintf('%s/output/movies/smoothing_error_in_image.png', floc), 'png');
ii=input('\nSave (tune vd, Rk ...) ? ([]=yes, other=no): ');
if(isempty(ii))
     fdatfile=[pathname, '/', name, '_F', '.mat'];
     save(fdatfile, 'Xh', 'P', 'Xi');
else
    fdatfile=[];
end
close(1);
