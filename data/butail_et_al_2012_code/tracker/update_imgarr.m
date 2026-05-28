function imgarr=update_imgarr(imgarr, k, k0, noise_std, flist, frmloc, optTrack)
%function imgarr=update_imgarr(imgarr, k, optTrack.br0, noise_std, flist, frmloc)
%
% function updates the imagearr incrementally by reading off images from
% the list of files

% if this is the first image then fill up the imgarr
if (k==k0)
    % fill up all that is possible in the begining and everything in
    % the end
    jj1=1;
    for jj=k-optTrack.br0:k+optTrack.br0
        img=imread([frmloc, flist(jj).name]);
        if(size(img,3)>1), img=rgb2gray(img); end
        if(noise_std>0), img=filter2(fspecial('gaussian', [3 3], noise_std), img); end
        % index is such that we want k to be right on optTrack.br0+1,
        % therefore we go back optTrack.br0 and add jj and then -1
        imgarr(:,:,jj1)=img;
        jj1=jj1+1;
    end
else
    % do a circshift 
    imgarr=circshift(imgarr, [0,0,-1]);

    % populate the last index with k+optTrack.br0 image
    img=imread([frmloc, flist(k+optTrack.br0).name]);
    if(size(img,3)>1), img=rgb2gray(img); end
    if(noise_std>0), img=filter2(fspecial('gaussian', [3 3], noise_std), img); end        
    imgarr(:,:,2*optTrack.br0+1)=img;
end
