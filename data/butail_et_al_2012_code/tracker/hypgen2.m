function hypmat = hypgen2(vmat, tids)
% function hypmat = hypgen(vmat)
% vmat is the validation matrix with t+1 columns and nz rows. 
% NOTE: the first column should always be 1s only denoting false alarms
%
% this function generates the hypotheses for MHT
% nz is the number of measurements at the current time step
%
% ref: ﻿[1] D. Reid, “An Algorithm for Tracking Multiple Targets,” IEEE
% Trans. on Automatic Control, 1979.
%
% Examples:
% 
% hypgen([1 1 1])
% hypgen([1 1 1; 1 0 1])
% hypgen([1 1 1; 1 0 1; 1 0 1]) gives 28 hypotheses
%
% tic; cs=5; hypgen([ones(cs,1), rand(cs)>.5]); toc % to test perf.
%
% Mar 2011, Sachit Butail

global gf;
global idl;

nz=size(vmat,1);

if nz > 10
    error('[!] Too many measurements. Have you performed clustering..?');
end

% add the ones in the first column denoting all measurements can be false
% vmat=[ones(nz,1) vmat];
eid=size(vmat,2)-1;

% new targets
vmat=[vmat eye(nz)];

% first measurment
hyp_jj=find(vmat(1,:)~=0)';

% second measurement onwards
for jj =2:nz
    % find the associated targets
    hh = find(vmat(jj,:)~=0)';
    % append into the hypotheses from previous measurement just as in Fig.
    % 2 of the paper
    hyp_jj=[kron(ones(size(hh,1),1), hyp_jj), kron(hh, ones(size(hyp_jj,1),1))];
end
        
% removing invalid hypotheses entries--meaning one target sends two measurements
nh=size(hyp_jj,1);

% just a flag vector to find valid hypotheses
val_hyp=zeros(nh,1);

% cycle through all hypotheses
for jj=1:nh
    % number of false alarms
    nf=sum(hyp_jj(jj,:)==1);
    
    % don't double count the false alarm
    if nf > 0, nf=nf-1; end
    
    if numel(unique(hyp_jj(jj,:)))+nf == nz
        val_hyp(jj)=1;
    end
end
    
% to be consistent with the paper terminology that 0 means a false alarm.
% if you want to be consistent with JPDA style then add 1 to all after you
% output
hypmat=hyp_jj(val_hyp~=0,:)-1;

hypmat1=hypmat;
for jj=1:eid
    hypmat(hypmat1==jj)=tids(jj);
end
