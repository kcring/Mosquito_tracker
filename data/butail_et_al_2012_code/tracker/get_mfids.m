function [male_id female_id]=get_mfids(floc)

metalist=dir([floc, 'calib/metadata_*.txt']);
if(size(metalist,1))
    metadata=csvread([floc, 'calib/', metalist(1).name]);
    if numel(metadata) ==5
        male_id=metadata(5);
        female_id=metadata(4);
    else
        fprintf('[I] metadata does not have male/female id appended. use append_metadata_with_mqids.m\n');
        male_id=0;
        female_id=0;
        return;
    end
end