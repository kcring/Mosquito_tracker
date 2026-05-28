function ret=fwback_diff(img_info, k, ball_radius, filter_noise_dev, max_min, file_list, file_loc, ud, ud_cam)
%function fg=fwback_diff(img_info, k, ball_radius, binary_t, area_t,
% filter_noise_dev, max_min, file_list, file_loc)
% ud =1 if you wish to use undistorted images
% ud_cam is cam structure that has to be passed if ud=1
% ret is a structure with ret.fg as the foreground and ret.bg as the
% background

image_width=img_info.Width;
image_height=img_info.Height;


bitval=2^img_info.BitDepth-1;

imgarr=zeros(image_height,image_width, 2*ball_radius+1);

jj=1;
for kk=k-ball_radius:k+ball_radius

    % get the current image
    Xi0=imread(strcat(file_loc,file_list(kk).name));

    if(size(Xi0,3)>1), Xi=rgb2gray(Xi0); else Xi=Xi0; end
    
    if(ud)
        Xi=undistort_image(Xi, ud_cam);
    end

    if(filter_noise_dev>0), Xi=filter2(fspecial('gaussian', [3 3], filter_noise_dev), Xi); end

    % populate the array
    imgarr(:,:,jj)=Xi;

    jj=jj+1;
end

% get the background
% depending on how you see the targets this could be max or min
if(max_min==1)
    bg=max(imgarr,[],3);
else
    bg=min(imgarr,[],3);
end
fg=imsubtract(bg, double(imgarr(:,:,ball_radius+1)));
fg=fg/bitval;
fg_raw=fg;

ret.bg=bg;
% this part can be used to vary the threshold on the fly
ret.fg_raw=fg_raw;