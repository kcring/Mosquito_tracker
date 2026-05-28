function mstr = extract_M(floc, camid_list)

for ii =1:size(camid_list,1)
    load([floc, '/output/data/', 'fg_', camid_list(ii,:), '.mat']);
    mst(ii).frm=mstr.frm;
    mst(ii).file=mstr.file;
    clear mstr
end
mstr=mst;