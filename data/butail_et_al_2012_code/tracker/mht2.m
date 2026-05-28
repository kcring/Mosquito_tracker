function mht2

% ---------------------------------- parameters, verify this!!
global params k tid

k=1; 

tid=1;

params.t_win=7; % sliding window
params.sigma_ep=diag([2 2]); % end points
params.sigma_w=100; % m^2/s^4
params.dt=1/10;
params.t_gate=25; 
params.t_area=[20, 150]; % min max blob areas
params.max_per_cluster=100; % a large value will suppress clustering for now..
params.scanback=2; % scanback
params.Np=20; % number of particles
params.P_D=0.9; % probability of detection
params.beta_f=1/100; % density of false targets
params.beta_n=1/100; % density of new targets
params.hyp_sort=3; % number of hypotheses to keep based on probability

% -------------------------------- initialize
Z{k}=[];
clustr=init_clustr;
tic
while k < 20
    
    fprintf('%d(%.1f s)..', k, toc);
    tic
    if ~mod(k,5), fprintf('\n'); end
    % Extract measurements from images (image processing)
%     blobs{k,1}=extract_measurements();

    % Using predicted positions of *confirmed* targets within clusters to
    % look for missing measurements by adaptively thresholding the region
    % in the image where you expect to find it
%     blobs{k,1}=find_missing_measurements();

    % Use epipolar constraint to convert blobs in the image pairs to 3D
    % measurements
    % For dry tests, pass preset measurements here..
    Z{k}=validate();

    % A cluster consists of hypotheses. A hypotheses is a measurement target
    % assignment. As time progresses, depending on
    % the scanback, the number of hypotheses at a given step is left to one
    % or more. 
    % The cluster looks like 
    % clustr(1)=struct(zidx, tidx, hypotheses(k), Xp), where hypotheses is
    % a further structure with fields 'value', 'costmat', 'prob'
    %
    % Note that the number of hypotheses in the first time step is empty.
    % cluster function also combines and separates clusters and the 
    % corresponding hypotheses. For clusters that are combined, the 
    % resulting hypotheses are a product of the number of hypotheses 
    % in each cluster.
    clustr=cluster(Z{k}, clustr);

    for ii_cl=1:size(clustr,2)
        % Hypotheses are computed for each cluster individually, as if that is
        % the only set of targets and measurements. In the first step, since
        % there are no targets, nothing happens here
        clustr(ii_cl)=compute_hypotheses(clustr(ii_cl), Z{k}); 

        % Hypothesis reduction is performed by N-scanback and thresholding. For
        % N-scanback, use the probabilities at the current timestep of all
        % hypotheses to select the one that is the most probable in the previous time step.
        % Then remove all the other parents and their child hypotheses. From the ones
        % remaining, sort according to probability and further reduce to a few
        % high ones. In the first step, nothing to be done since there are no
        % targets.
        clustr(ii_cl)=hypotheses_reduction(clustr(ii_cl)); 
        
        % For each unassigned measurement, select a target. In the first step,
        % this is already done for each cluster
        clustr(ii_cl)=initialize_and_update(Z{k}, clustr(ii_cl));

        clustr(ii_cl)=predict(clustr(ii_cl));
        
    end
    
    % update the time step
    k=k+1;
end

% plot and compare
figure(1); gcf; clf;
z=cat(2, Z{:});
plot(z(1,:), z(2,:), 'kx'); 
hold on;
for ii_cl=1:size(clustr,2)
    ids=unique(clustr(ii_cl).Xp(:,2))';
    ids=ids(ids~=0);
    clrs=colormap(jet(numel(ids)));
    for tgt_ii=ids
        track=zeros(3,k-1);
        for k_ii=1:k-1
            track(1:2,k_ii)=postest(clustr(ii_cl).Xp(clustr(ii_cl).Xp(:,2)==tgt_ii & ...
                                clustr(ii_cl).Xp(:,1)==k_ii, 4:5), [], 1);
        end
        plot(track(1,:), track(2,:), 'o-', 'color', clrs(find(tgt_ii==ids),:));
    end
end
keyboard
% ----------------------------------
function blobs=extract_measurements(frmloc, identifier1, identifier2)

global params k

global imgarr

bitval=8; % usually but some images are more


% initialize
for cc=1:numel(identifier2)
    imstream(cc).flist=dir([frmloc, identifier1, identifier2,'*.*']);
    bgparams=init_bgparams(imstream(cc), cc, floc, frmloc, ...
                    expected_nmq, optTrack);
end


function find_missing_measurements()

function zk=validate()
global params k

load('./vicsek2dz_N5.mat');

% low noise
nz=size(x1,1);
zk=[x1(:,k)'; x2(:,k)']+randn(2,nz)*.01;

% high noise
% nz=size(x1,1);
% zk=[x1(:,k)'; x2(:,k)']+randn(2,nz)*.1;


% skipping measurements
ns=rand;
if ns < 0 % 0.25
    ns=ceil(ns*nz);
    rm_z=randperm(nz);
    rm_z=rm_z(1:ns);
    zk(:,rm_z)=[];
end
% missing measurements

function clustr=cluster(Z, clustr)
% NOTE: clusters are target-driven in contrast to the general MHT algorithm
% which is measurement driven
% Use gating volume of each target within a cluster to add measurements to 
% that cluster. A cluster is the smallest set of measurements and targets 
% that exist independently; combine/divide existing clusters as needed.
%
% this function takes all the measurements at current timestep, clusters at
% previous timestep, and hypotheses at previous time step to create new
% clusters. if no clusters exist at current time step, then do k-means such
% that there are on an average 5 targets per cluster

global params k

% In the first time-step simply use k-means and set one-to-one measurement
% target assignment. 
if k ==1
    nz=size(Z,2);
    if round(nz/params.max_per_cluster) > 1
        idx=kmeans(Z', round(nz/params.max_per_cluster));
    else
        idx=1;
    end
    for ii_cl=unique(idx)'
        
        clustr(ii_cl)=init_clustr;
        clustr(ii_cl).zidx{k}=find(idx==ii_cl)';
    end
else
    nz=size(Z,2);
    ncl=size(clustr,2);
    % construct a validation matrix [nzxnc] for clusters and measurements
    % then
    % construct a validation matrix for each cluster i.e. use target gating
    % volumes to find which measurements are to be assigned to that
    % cluster. Note this is different from k-means but the output should be
    % the same, that is the clustr.zidx will have the measurement numbers
    % of that time-step. After that it is all the same.
    
    validation_matrix_Z_clustr=zeros(nz, ncl);
    
    for ii_cl=1:size(clustr,2)
        
        % initialize
        zidx=[]; 
        
        % construct the covariance for each projected estimate , which is
        % in each hypothesis
        nodes_last_step=find(clustr(ii_cl).node ==k-1);
        for hyp_ii=1:numel(nodes_last_step) 
            node_id=nodes_last_step(hyp_ii);
            
            hypothesis=clustr(ii_cl).hypotheses.get(node_id);
            icov=zeros(nz, numel(hypothesis));
            for tgt_ii=1:numel(hypothesis)
                target_id=hypothesis(tgt_ii);
                rows2=clustr(ii_cl).Xp(:,3)==node_id & ...
                      clustr(ii_cl).Xp(:,2)==target_id & ...
                      clustr(ii_cl).Xp(:,1)==k;
                
                % covariance *of the measurement estimate*
                S=cov(clustr(ii_cl).Xp(rows2, 4:6));
                
                zh=postest(clustr(ii_cl).Xp(rows2, 4:6), [], 1);
                
                for z_ii=1:nz
                    % *TODO* check the 2D
                    icov(z_ii, tgt_ii)=...
                        gatefun(zh(1:2)', Z(:,z_ii), S(1:2,1:2));
                end
                zidx=[zidx, find(icov(:,tgt_ii)<params.t_gate)'];
                
            end
            
            % two hypotheses may point to the same measurement and thus
            % multicount the measurement. 
            clustr(ii_cl).zidx{k}=unique(zidx);
            
            % now create the costmatrix based on number of measurements and
            % targets in this hypothesis
            clustr(ii_cl).costmat=clustr(ii_cl).costmat.set(node_id, icov(clustr(ii_cl).zidx{k}, :));

        end
        
        validation_matrix_Z_clustr(clustr(ii_cl).zidx{k}, ii_cl)=1;
    end
    
    % is there a measurement that does not belong to any cluster--then make
    % a new cluster
    zidx_no_clustr=find(sum(validation_matrix_Z_clustr,2)==0,1);
    if ~isempty(zidx_no_clustr)
        for z_ii=zidx_no_clustr
            ii_cl=ii_cl+1;
            clustr(ii_cl)=init_clustr;
            clustr(ii_cl).zidx{k}=z_ii;
        end
    end
        
    % is a measurement in two or more clusters--combine them
%     zidx_multi_clustr=find(sum(validation_matrix_Z_clustr,2)>1,1);
%     if ~isempty(zidx_multi_clustr)
%         keyboard
%         ii_cl=ii_cl-1;
%     end
    
    % is there a target that is assigned to a single measurement within a
    % cluster. for e.g. is the validation matrix of the cluster such that
    % there is a target that is only assigned to measurement 3 and false
    % target, then split that cluster
    
end


function cl=init_clustr()

global params

% Xp is Np x targets x hypotheses
% hypotheses is a tree structure
% nodes are updated so that each timestep is a level
% costmat corresponds to each hypothesis 

cl=struct('zidx', [], ...
            'Xp', zeros(params.Np*10*10, 1+1+1+3+3), ...
            'hypotheses', tree, 'node', tree, 'costmat', ...
            tree, 'prob', tree);


function clustr_ii=compute_hypotheses(clustr_ii, Zk)

global params k tid

% if this is a new clustr
if isempty(clustr_ii.node.Node{1})
    nz_cl=numel(clustr_ii.zidx{k});
    targets=(tid:tid+nz_cl-1);
        
    % hypotheses are just the targets in that order
    clustr_ii.hypotheses=tree(targets);
    
    clustr_ii.node=tree(k);
    clustr_ii.prob=tree(1);
    clustr_ii.costmat=tree(eye(numel(targets)));

    % update the target ids
    tid=tid+nz_cl;
else
    nodes_last_step=find(clustr_ii.node ==k-1);
    
    % check that all trees are in sync
    if ~clustr_ii.hypotheses.issync(clustr_ii.costmat) || ...
       ~clustr_ii.hypotheses.issync(clustr_ii.prob) || ...
       ~clustr_ii.hypotheses.issync(clustr_ii.node)
        error('[%d] cluster trees not in sync.. (compute_hypotheses)', k);
    end
        
    % for every hypothesis in the last time step
    for hyp_ii=1:numel(nodes_last_step) 
        parent_node=nodes_last_step(hyp_ii);
        hypothesis=clustr_ii.hypotheses.get(parent_node);
        costmat=clustr_ii.costmat.get(parent_node);
        vmat=[ones(size(costmat,1),1), costmat<params.t_gate];
        try
            new_hypotheses=hypgen2(vmat, hypothesis, costmat);
        catch
            keyboard
        end
        % probability sum to noramlize
        psum=0;
        for new_hyp_ii=1:size(new_hypotheses,1) 
           [clustr_ii.hypotheses, node1]=clustr_ii.hypotheses.addnode(parent_node, ...
                        new_hypotheses(new_hyp_ii, :));
            clustr_ii.node=clustr_ii.node.addnode(parent_node, k);
            clustr_ii.prob=clustr_ii.prob.addnode(parent_node, 0);
            clustr_ii.costmat=clustr_ii.costmat.addnode(parent_node, []);
            parent.node=parent_node;
            parent.assignment=clustr_ii.hypotheses.get(parent_node);
            parent.prob=clustr_ii.prob.get(parent_node);

            clustr_ii.prob=clustr_ii.prob.set(node1,hypothesis_prob(new_hypotheses(new_hyp_ii,:), ...
                        parent, clustr_ii.Xp, Zk(:, clustr_ii.zidx{k})));
            psum=psum+clustr_ii.prob.get(node1);    
        end
        
        % normalize the probability
        for node_ii=clustr_ii.node.getchildren(parent_node)
            prob_ii=clustr_ii.prob.get(node_ii);
            clustr_ii.prob=clustr_ii.prob.set(node_ii, prob_ii/psum);
        end
    end
end

function hc=hypothesis_prob(hca, hp, Xp, Zk)

global k params

hpa=hp.assignment;

%target indices of those targets that are existing
t_d=hca(ismember(hca, hpa));
t_d=t_d(t_d~=0);

% Number of measurements associated with existing targets 
N_d=numel(t_d);


% Number of measurements associated false targets
N_f=sum(hca==0);

% Number of measurements associated with new targets
N_n=numel(hca)-N_d-N_f;

% number of existing targets from previous hyp
N_t=sum(hpa>0);

% Equation 16 in paper computing probability of existing targets
prob_detection=1;
for t=1:N_d
    target_id=t_d(t);

    rows2=Xp(:,3)==hp.node & ...
          Xp(:,2)==target_id & ...
          Xp(:,1)==k-1;
 
    % covariance *of the measurement estimate*

    S=cov(Xp(rows2, 4:6));
    
    zh=postest(Xp(rows2, 4:6), [], 1);
    
    % hca==tid is the measurement id 
    % *TODO* 2D to 3D
    prob_detection=prob_detection*normal2(Zk(:,hca==target_id), zh(1:2)', S(1:2,1:2));
end

% equation 16 in Reid's paper
hc=params.P_D^N_d*(1-params.P_D)^(N_t-N_d)*params.beta_f^N_f*...
                params.beta_n^N_n*prob_detection*hp.prob;

function val=normal2(x,mu,covar)
% x    : d x n; n data points with dim d
% mu   : d x 1
% covar: d x d covariance matrix
% Copyright 2009, Ali Bahramisharif
% This code is free to change, use and re-distribute.

[d, n]=size(x);
x = x-mu*ones(1,n);
x = x';

val = exp(-0.5*sum((x/(covar)).*x, 2))/ sqrt((2*pi)^d*abs((1e-10) +det(covar)));

function clustr_ii=hypotheses_reduction(clustr_ii)
global params k

% sort reduction, *TODO* make it less conservative
if k > params.scanback
    % get all parent hypotheses from previous time step
    nodes_last_step=find(clustr_ii.node ==k-1);
    
    for hyp_ii=1:numel(nodes_last_step) 
        parent_node=nodes_last_step(hyp_ii);
        
        % get all the children (these are at the current time-step)
        children=clustr_ii.node.getchildren(parent_node);
        prob_cc=zeros(1,numel(children));
        for node_cc_ii=1:numel(children)
            prob_cc(node_cc_ii)=clustr_ii.prob.get(children(node_cc_ii));
        end
        
        [~, idx]=sort(prob_cc);
        n_nodes=numel(idx);
        nodes_to_remove=children(idx(1:n_nodes-params.hyp_sort));
        
        % sort these as well since we remove from the end of the tree in a loop
        nodes_to_remove=sort(nodes_to_remove, 'descend');
        if n_nodes > params.hyp_sort
            for node_ii=nodes_to_remove
                clustr_ii.hypotheses=clustr_ii.hypotheses.removenode(node_ii);
                clustr_ii.prob=clustr_ii.prob.removenode(node_ii);
                clustr_ii.node=clustr_ii.node.removenode(node_ii);
                clustr_ii.costmat=clustr_ii.costmat.removenode(node_ii);
            end


            % normalize the probability
            psum=0;
            for node_ii=clustr_ii.node.getchildren(parent_node)
                psum=clustr_ii.prob.get(node_ii)+psum;
            end
            for node_ii=clustr_ii.node.getchildren(parent_node)
                prob_ii=clustr_ii.prob.get(node_ii);
                clustr_ii.prob=clustr_ii.prob.set(node_ii, prob_ii/psum);
            end
        end
    end
    
    % check that all trees are in sync
    if ~clustr_ii.hypotheses.issync(clustr_ii.costmat) || ...
       ~clustr_ii.hypotheses.issync(clustr_ii.prob) || ...
       ~clustr_ii.hypotheses.issync(clustr_ii.node)
        error('[%d] cluster trees not in sync.. (hypotheses_reduction)', k);
    end
end


function icov=gatefun(x,z,S) 

icov=(x-z)'/S*(x-z);


function hypmat = hypgen2(vmat, tids, cost)
% function hypmat = hypgen(vmat)
% vmat is the validation matrix with t+1 columns and nz rows. 
% NOTE: the first column should always be 1s only denoting false alarms
%
% this function generates the hypotheses for MHT
% nz is the number of measurements at the current time step
%
% ref: [1] D. Reid, An Algorithm for Tracking Multiple Targets,” IEEE
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


global params tid

nz=size(vmat,1);

if nz > params.max_per_cluster
    % run global assignment algorithm if the measurements are too many
    % because then the MHT hypotheses will explode!
    % munkres gives the assignments based on measurements (in columns) to
    % targets 
    hypmat1=munkres(cost');
    
    [~, assign]=max(hypmat1, [], 1);
    
    assign(sum(hypmat1,1)==0)=0;

    hypmat=zeros(1, size(cost,1));
    for jj=1:numel(tids)
        hypmat(assign==jj)=tids(jj);
    end

elseif nz <= params.max_per_cluster && nz > 0

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
    try
    for jj=1:eid
        hypmat(hypmat1==jj)=tids(jj);
    end
    catch
        keyboard
    end
    % new ids
    ids_all=unique(hypmat1(:));
    new_ids=ids_all(ids_all>max(eid));
    for jj=new_ids'
        tid=tid+1;
        hypmat(hypmat1==jj)=tid;
    end
elseif nz==0
    hypmat=[];
end

function clustr_ii=initialize_and_update(Z, clustr_ii)

global params k 


% dummy initialization and update, fix it later **TODO**
kidx=find(clustr_ii.Xp(:,1)~=0);
if isempty(kidx), kidx=0; end
iter=kidx(end)+1;
% for each hypothesis
nodes_k=find(clustr_ii.node ==k);

for hyp_ii=1:numel(nodes_k)
    node_id=nodes_k(hyp_ii);
    
    % and each measurement assignment
    hypothesis=clustr_ii.hypotheses.get(node_id);
    parent_node=clustr_ii.hypotheses.Parent(node_id);

    for z_ii=1:numel(hypothesis)
        
        % new target or not
        rows=clustr_ii.Xp(:,2)==hypothesis(z_ii) & clustr_ii.Xp(:,1)==k ...
                & clustr_ii.Xp(:,3)==parent_node;
        if ~sum(rows)
        
            % note that now we will have a target in multiple hypotheses
            % according to the measurement assignment!
            % [k, tid, hyp#, ...]
            clustr_ii.Xp(iter:iter+params.Np-1,:)=ones(params.Np,1)*[k, ...
                            hypothesis(z_ii), node_id, ...
                            Z(1:2,clustr_ii.zidx{k}(z_ii))', 0, ...
                            1, 1, 0];
            % make some noise *TODO* change
            clustr_ii.Xp(iter:iter+params.Np-1,4:5)= ...
                clustr_ii.Xp(iter:iter+params.Np-1,4:5)+...
                randn(params.Np,2)*.05;
            
            iter=iter+params.Np;
        else
            clustr_ii.Xp(rows,:)=[];
            % *TODO* this should be different
            clustr_ii.Xp(iter:iter+params.Np-1,:)=ones(params.Np,1)*[k, ...
                            hypothesis(z_ii), node_id, ...
                            Z(1:2,clustr_ii.zidx{k}(z_ii))', 0, ...
                            1, 1, 0];
            % make some noise *TODO* change
            clustr_ii.Xp(iter:iter+params.Np-1,4:5)= ...
                clustr_ii.Xp(iter:iter+params.Np-1,4:5)+...
                randn(params.Np,2)*.05;
            iter=iter+params.Np;
        end
                  
    end
end

function clustr_ii=predict(clustr_ii)

global params k 


% again dummy for now but should work... **TODO**
kidx=find(clustr_ii.Xp(:,1)~=0);
if isempty(kidx), kidx=0; end
iter=kidx(end)+1;

% for each hypothesis
nodes_k=find(clustr_ii.node ==k);

for hyp_ii=1:numel(nodes_k)
    node_id=nodes_k(hyp_ii);
    rows=clustr_ii.Xp(:,3)==node_id;
    
    % this is where the motion model goes, for now we just say that it is
    % constant velocity ** NOTE that we are updating multiple targets as
    % they belong to this hypothesis
    clustr_ii.Xp(iter:iter+sum(rows)-1, 1:6)=[(k+1)*ones(sum(rows),1), ...
                                        clustr_ii.Xp(rows,2:3), ...
                                        clustr_ii.Xp(rows, 4:6)+clustr_ii.Xp(rows, 7:9)*params.dt];
    clustr_ii.Xp(iter:iter+sum(rows)-1, 7:9)=clustr_ii.Xp(rows, 7:9)+randn(sum(rows),3)*.05;
    
    % **TODO** delete this later
    clustr_ii.Xp(iter:iter+sum(rows)-1, 6)=0; clustr_ii.Xp(iter:iter+sum(rows)-1,9)=0;
    
    iter=iter+sum(rows);
end


function estx = postest(x, wts, flag)

%{
This function will return an estimate of a distribution, based on what flag is passed.
So, for e.g. if flag passed is 1 and x is a matrix, the function will compute the mean along
the longer dimension
if flag is 2 then the function will compute the max along the longer dimension... and so on.

ok. so there's a check now that the longer dimension has to be columns!!
%}

switch flag
    case 1
        estx = mean(x,1);
    case 2
        [~, idx] = max(wts);
        estx = x(idx,:);
    case 3
        estx = mode(x,1);
    otherwise
        error('brrrr!!');
end
    
