function reformat_data(datafile)

load(datafile)

if (exist('Xi', 'var') && Xi.nX ~=6) || (exist('nx', 'var') && nx~=6)
    
    if exist('Xi', 'var'), nx=Xi.nX; end
    if exist('nx', 'var'), nx=nx; end
    
    fprintf('[I] Reformatting data according to new setup....\n');
    nt=size(Xh,1)/nx;

    for ii=1:nt
        r0=getind(nx,1, ii, 1:6, 1);
        r1=getind(6,1,ii,1:6,1);
        Xh1(r1,:)=Xh(r0,:);
    end
    Xi=strXi;
    Xh=Xh1;
    
    if exist('Z', 'var')
        save(datafile, 'Xh', 'Xi', 'Z');
    else
        save(datafile, 'Xh', 'Xi');
    end
end