% join tracks

clear all;
% --- glue
% load('/media/mali2010rep/2010/Aug/21/NIH/tiff/20100821_NIH_19.02.08/output/data/data_mq_EXP_2010_08_21_NIH_jose.mat', 'Xh');
% data_jose.Xh=Xh;
% clear('Xh');
% load('/media/mali2010rep/2010/Aug/21/NIH/tiff/20100821_NIH_19.02.08/output/data/data_mq_EXP_2010_08_21_NIH_19.02.08.mat');
% Xi=strXi;
% show_tracked_mq(data_jose.Xh, Xi, 0, 'h', 1);
% fprintf('=============================================================\n');
% show_tracked_mq(Xh, Xi, 0, 'h', 1);
% id(1)=63; t(1)=1;
% id(2)=63; t(2)=759;
% 
% 
% 
% [r1 c1]=getind(Xi.nX, t(1), id(1), 1:Xi.nX, 1);
% [r2 c2]=getind(Xi.nX, t(2), id(2), 1:Xi.nX, 1);
% 
% Xh(r1,c2:end)=data_jose.Xh(r2,c2:end);
% 
% % Xh(r2,:)=0;
% clear('data_jose');
% save('/media/mali2010rep/2010/Aug/21/NIH/tiff/20100821_NIH_19.02.08/output/data/data_mq_EXP_2010_08_21_NIH_19.02.08.mat');


% --- swap tracks
load('/home/sachit/Dropbox/CDCL mosquito data/data/MalesOnly/21aug2010_19.02.08.707/output/data/data_mq_EXP_2010_08_21_NIH_19.02.08.mat');
id(1)=27;
id(2)=63;

tswap=1;
[r1 c1]=getind(Xi.nX, tswap, id(1), 1:Xi.nX, 1);
[r2 c2]=getind(Xi.nX, tswap, id(2), 1:Xi.nX, 1);

tXh=Xh(r2,c2:end);
Xh(r2,c2:end)=Xh(r1,c1:end);
Xh(r1,c1:end)=tXh;
save('/home/sachit/Dropbox/CDCL mosquito data/data/MalesOnly/21aug2010_19.02.08.707/output/data/data_mq_EXP_2010_08_21_NIH_19.02.08.mat');
