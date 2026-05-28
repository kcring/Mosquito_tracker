function [expname image_id]=scan_expfile(floc)


expfile=[floc, '/calib/', 'expfile.txt'];
if exist(expfile, 'file')
    fid=fopen(expfile);
    txt=textscan(fid, '%s%s');
    expname=char(txt{1});
    image_id=char(txt{2});
    if isempty(image_id)
        fprintf('[!] Could not find image_id in the expfile. Edit the file or run setupdir again\n');
        return;
    end
else
    fprintf('[!] Run setupdir first ...\n');
    expname=[];
    image_id=[];
    return;
end