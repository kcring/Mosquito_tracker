function append_metadata_with_mqids

% this function will append female, and focal male ids to the metadata
% file.

% select mosquito id
close all

addpath core/
run initproc

fprintf('I will now append the metadata file with important mqids... \nYou can rerun this function if you change the data...\n');

%% convert to world coordinates
fprintf('Loading metadata from %s ....\n', floc);



metafile=dir([floc 'calib/' 'metadata_*.txt']);
if(size(metafile,1))
    meta_file=[floc, 'calib/', metafile(1).name];
else
    fprintf('[!]Could not find a metadata file\n');
end

fprintf('Loading data file from this dataset...\n');
load([floc, 'output/data/', sprintf('data_mq_%s.mat', expname)], 'Xh');

%% hand tracked mosquitoes.... select mqB, mqA
m_flags=Xh(nf(2):nx:end,:);
m_ind=find(m_flags==2);
[m_mqid k]=ind2sub(size(m_flags), m_ind);

m_mqid=unique(m_mqid);
if(~isempty(m_mqid))
    fprintf('Total trackone tracks = %d ... \n', numel(m_mqid));
    fprintf('\n');    
    fprintf('id \t start \t end\n');
    fprintf('---------------------\n');
    for mq_ii=m_mqid'
        [r c]=getind(nx, 1, mq_ii, 1, 1);
        nz_ind=find(Xh(r,:)~=0);
        fprintf('%d \t %d \t %d\n', mq_ii, nz_ind(1), nz_ind(end));
    end
    fprintf('=====================\n');    
end

meta_info=csvread(meta_file);

if numel(meta_info) > 3
    curr_female_mq_id=meta_info(4);
    curr_focal_male_mq_id=meta_info(5);
else
    curr_female_mq_id=0;
    curr_focal_male_mq_id=0;
end
female_id=input(sprintf('Current female mq id = %d. []= keep, newid=newid: ', curr_female_mq_id));
if ~isempty(female_id)
    meta_info(4)=female_id;
end
focal_male_id=input(sprintf('Current focal male mq id = %d. []= keep, newid=newid: ', curr_focal_male_mq_id));
if ~isempty(focal_male_id)
    meta_info(5)=focal_male_id;
end

csvwrite(meta_file, meta_info);
fprintf('done....\n');

    
  