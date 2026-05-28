function setupdir

run config

exp_dir=floc;
if ~exist(exp_dir, 'dir')
    fprintf('Did not find %s ..\n', exp_dir);
    resp=input('Do you want to create it (y/n)?', 's');
    if strcmp(resp, 'y')
        mkdir(exp_dir);
    else
        return;
    end
end

mkdir(exp_dir, '/frames');
mkdir(exp_dir, '/output');
mkdir([exp_dir, '/output'], '/data');
mkdir([exp_dir, '/output'], '/movies');
mkdir(exp_dir, '/calib');

% check and fix folders and set status as 'x' if done
% display 

% metadata
% number of fish, tank dimensions [ L X W X H] 
metafile_list=dir([floc, '/calib/metadata*.txt']);
if size(metafile_list,1)
    metafile=metafile_list(1).name;
    fprintf('[I] %s file exists... re-editing file\n', metafile);
    md=csvread([floc, '/calib/', metafile]);
else
    md=[0,0,500];
end
kk=input(sprintf('Azimuth / Compass direction (deg.) []=%d: ',md(1)));
if ~isempty(kk), md(1)=kk; end

kk=input(sprintf('Inclination (deg.) []=%d: ', md(2)));
if ~isempty(kk), md(2)=kk; end

if numel(md)==2, md(3)=500; end % backwards compatibility
kk=input(sprintf('Camera height (mm) []=%d: ', md(3)));
if ~isempty(kk), md(3)=kk; end

% naming etc.
expfile=[exp_dir, '/calib/', 'expfile.txt'];
if exist(expfile, 'file')
    fprintf('[I] expfile file exists...Delete it if you want to reset the values \n');
    [expname image_id]=scan_expfile(exp_dir);
else
    [frmloc image_id exp_name]=getImageId;
    fid=fopen(expfile, 'w');
    fprintf(fid, '%s %s', exp_name, image_id);
    fclose(fid);
end


csvwrite([floc, '/calib/metadata_', image_id, '.txt'], md);