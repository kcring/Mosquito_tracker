function bgparams=init_bgparams(imstream, camid, floc, frmloc, expected_nmq, optTrack)
% function bgparams=init_bgparams(imstream, camid, floc, frmloc, expected_nmq, optTrack)
% function to automatically compute background parameters.. right now only
% binary_t is set more later

fprintf('[I] Computing (one time only) background parameters for camera %d...\n', camid);


bgparams.area_t=[10 300]; % min max 
bgparams.br=3;
bgparams.noise_std=1;
bgparams.expected_nmq=expected_nmq*2;
bgfile=sprintf('%scalib/cam%d_bgparams.csv', floc, camid);


if ~exist(bgfile, 'file')

    nframes_test=9;
    % nframes=size(imstream.flist);
    % 
    % step1=floor(nframes/100);
    k0=optTrack.br0+1;

    figure(1); gcf; clf;
    imshow([frmloc, imstream.flist(1).name]);
    xlabel('Select region of interest. Use this to remove trees and other swarms...', 'Color', 'r');
    roi=getrect;
    bgparams.roi=[ceil(roi(1)) ceil(roi(2)) floor(roi(3)) floor(roi(4))];
    binary_test=.03*ones(nframes_test+1,1);
    
    jj=1;
    for k=k0:k0+nframes_test

        fprintf('.');

        % update image array
        imstream.imgarr=update_imgarr(imstream.imgarr, k, k0, ...
                        bgparams.noise_std, imstream.flist, frmloc, optTrack);
        Zk=[];
        
        while size(Zk,2) < expected_nmq
            imstream.area_t=bgparams.area_t;    
            imstream.binary_t=binary_test(jj);    
            imstream.br=bgparams.br;
            imstream.roi=bgparams.roi;
            
            Zk=getZ(imstream, -1, optTrack);
            binary_test(jj)=binary_test(jj)*.97;
        end
        jj=jj+1;
    end

    bgparams.binary_t=mean(binary_test);
    fid=fopen(bgfile, 'w');
    fprintf(fid, '%.4f, %.1f, %d, %.1f, %d, %d, %d, %d, %d', bgparams.binary_t, bgparams.area_t(1), ...
                        bgparams.br, bgparams.noise_std, bgparams.expected_nmq, bgparams.roi(1), ...
                        bgparams.roi(2), bgparams.roi(3), bgparams.roi(4));
    fclose(fid);
else
    fprintf('[I] found background params file.. reading..');
    bgread=csvread(bgfile);
	bgparams.binary_t=bgread(1);
    bgparams.area_t(1)=bgread(2);
    bgparams.br=bgread(3);
    bgparams.noise_std=bgread(4);
    bgparams.expected_nmq=bgread(5);
    bgparams.roi=[bgread(6), bgread(7), bgread(8), bgread(9)];
end

fprintf('done. The value for binary_t=%.4f\n', bgparams.binary_t);
