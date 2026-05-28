function compareAutoWithManual

run config
Xi=strXi;

lw=2;
fs=16;
%% load files
if isunix
    homedir='/home/sachit/';
elseif ispc
    homedir='C:\Users\Sachit\';
end

data(1).mdir=[homedir 'Dropbox/CDCL mosquito data/data/MalesOnly/29aug2010_CDCL_19.00.07.536/'];
data(1).adir=[homedir 'SampleData/mosquitoes/20100829_CDCL_19.00.07/'];

data(1).name='29aug2010_CDCL_19.00.07.536';
data(1).nframes=750;

data(2).mdir=[homedir 'Dropbox/CDCL mosquito data/data/MalesOnly/26aug2010_NIH_18.58.42.588/'];
data(2).adir=[homedir 'SampleData/mosquitoes/20100826_NIH_18.58.42/'];

data(2).name='26aug2010_NIH_18.58.42.588';
data(2).nframes=750;



for ii=1:size(data,2)
    fprintf('Reading data for %s..\n', data(ii).mdir);
    metalist=dir([data(ii).mdir, 'calib/metadata_*.txt']);
    if(size(metalist,1))
         expfile=[data(ii).mdir, '/calib/', 'expfile.txt'];
        if exist(expfile, 'file')
           expname=scan_expfile(data(ii).mdir);
        else
            fprintf('[!] Did not find expfile for this dataset ...\n');
            return;
        end
        if exist([data(ii).mdir, 'output/data/', sprintf('data_mq_%s.mat', expname)], 'file');
            data(ii).datafile=[data(ii).mdir, 'output/data/', sprintf('data_mq_%s.mat', expname)];
            load(data(ii).datafile, 'Xh');
%             nz_ind=find(Xh(3,:)~=0);
            data(ii).mXh=Xh;%(any(Xh,2),1:nz_ind(end));
            clear Xh;
        end
    end
    
    fprintf('Reading data for %s..\n', data(ii).adir);
    metalist=dir([data(ii).adir, 'calib/metadata_*.txt']);
    if(size(metalist,1))
         expfile=[data(ii).adir, '/calib/', 'expfile.txt'];
        if exist(expfile, 'file')
           expname=scan_expfile(data(ii).adir);
        else
            fprintf('[!] Did not find expfile for this dataset ...\n');
            return;
        end
        if exist([data(ii).adir, 'output/data/', sprintf('data_mq_auto_%s.mat', expname)], 'file');
            data(ii).datafile=[data(ii).adir, 'output/data/', sprintf('data_mq_auto_%s.mat', expname)];
            load(data(ii).datafile, 'Xh');

            data(ii).aXh=Xh;%(any(Xh,2),1:nz_ind(end));
            clear Xh;
        end
    end
end

% compute mean/std
for ii=1:size(data,2)
     mr(1).v=data(ii).mXh(1:Xi.nX:end,:);
     mr(2).v=data(ii).mXh(2:Xi.nX:end,:);
     mr(3).v=data(ii).mXh(3:Xi.nX:end,:);
     mnmq(ii).v=sum(mr(3).v~=0);
     
     ar(1).v=data(ii).aXh(1:Xi.nX:end,:);
     ar(2).v=data(ii).aXh(2:Xi.nX:end,:);
     ar(3).v=data(ii).aXh(3:Xi.nX:end,:);
     anmq(ii).v=sum(ar(3).v~=0);
     
     
     % means
     for jj=1:3
         data(ii).mr(jj).m=sum(mr(jj).v)./mnmq(ii).v;
         data(ii).ar(jj).m=sum(ar(jj).v)./anmq(ii).v;
     end
         
     % variance
     for jj=1:3
         mzi=find(mr(jj).v==0);
         azi=find(ar(jj).v==0);
         
         mdiff=(mr(jj).v-ones(size(mr(jj).v,1),1)*data(ii).mr(jj).m);
         mdiff=mdiff.^2; mdiff(mzi)=0;
         data(ii).mr(jj).std=sqrt(sum(mdiff)./mnmq(ii).v);
         
         adiff=(ar(jj).v-ones(size(ar(jj).v,1),1)*data(ii).ar(jj).m);
         adiff=adiff.^2; adiff(azi)=0;
         data(ii).ar(jj).std=sqrt(sum(adiff)./anmq(ii).v);
     end
end


% centroid
show=0;
if show
for ii=1:size(data,2)
    figure(ii); gcf; clf;

    nf=data(ii).nframes;
    alim=[-5 5;
       -5 5;
       20 30;
       0 20];

    for jj=1:3
         subplot(4,1,jj); gca;
         % -- manual
         plot((1:nf)/25, data(ii).mr(jj).m(1:nf)/100, 'r', 'LineWidth', lw);
         hold on;

         % -- auto       
         plot((1:nf)/25, data(ii).ar(jj).m(1:nf)/100, 'b', 'LineWidth', lw);

         set(gca, 'fontsize', fs);   
         ylabel(sprintf('r_%d (cm)',jj));
         set(gca, 'ylim', alim(jj,:));
         box off
         if jj==1
             title([data(ii).name ' (Camera frame)'], 'Interpreter', 'none');
         end
    end 
    jj=4;
    subplot(4,1,jj); gca;
    set(gca, 'fontsize', fs);
    plot((1:nf)/25, mnmq(ii).v(1:nf), 'r', 'linewidth', lw);
    hold on;
    plot((1:nf)/25, anmq(ii).v(1:nf), 'b', 'linewidth', lw);
    ylabel('# of mq');
    set(gca, 'ylim', alim(jj,:));
    legend('manual', 'auto');
    box off
    xlabel('sec');


    if 1
     set(gcf,'PaperPositionMode','auto')
     print('-depsc', sprintf('/tmp/%s.eps', data(ii).name));
    end
end
end
 
show=1;
if show
for ii=1:size(data,2)
    nf=data(ii).nframes;
    alim=[-5 5;
       -5 5;
       23 33;
       0 20];
    figure(ii); gcf; clf;

    for jj=1:3
        subplot(4,1,jj); gca;

        boundary1=data(ii).mr(jj).m(1:nf)/100+data(ii).mr(jj).std(1:nf)/100;
        boundary2=data(ii).mr(jj).m(1:nf)/100-data(ii).mr(jj).std(1:nf)/100;
        ts=(1:nf)/25;  
        nani=find(~isnan(boundary1));
        fill([ts(nani), fliplr(ts(nani))], [boundary1(nani), fliplr(boundary2(nani))],...
           'r', 'FaceAlpha', .25, 'EdgeColor', [1 0 0]+.75*(ones(1,3)-[1 0 0]));

        hold on;   
        boundary1=data(ii).ar(jj).m(1:nf)/100+data(ii).ar(jj).std(1:nf)/100;
        boundary2=data(ii).ar(jj).m(1:nf)/100-data(ii).ar(jj).std(1:nf)/100;
        nani=find(~isnan(boundary1));
        fill([ts(nani), fliplr(ts(nani))], [boundary1(nani), fliplr(boundary2(nani))],...
           'b', 'FaceAlpha', .25, 'EdgeColor', [0 0 1]+.75*(ones(1,3)-[0 0 1]));
        set(gca, 'fontsize', fs);   
        ylabel(sprintf('r_%d (cm)',jj));
        set(gca, 'ylim', alim(jj,:));
        box off
        if jj==1
            title(sprintf('standard deviation \n%s (Camera frame)', data(ii).name), 'Interpreter', 'none');
        end       
    end
    
    jj=4;
    subplot(4,1,jj); gca;
    set(gca, 'fontsize', fs);
    plot((1:nf)/25, mnmq(ii).v(1:nf), 'r', 'linewidth', lw);
    hold on;
    plot((1:nf)/25, anmq(ii).v(1:nf), 'b', 'linewidth', lw);
    ylabel('# of mq');
    set(gca, 'ylim', alim(jj,:));
    legend('manual', 'auto');
    box off
    xlabel('sec');
    if 0
        set(gcf,'PaperPositionMode','auto')
        print('-depsc', sprintf('/tmp/%s.eps', data(ii).name));
    end
end
end