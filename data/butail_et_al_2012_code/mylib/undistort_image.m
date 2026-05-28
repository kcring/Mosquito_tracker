function iundist = undistort_image(im, cam, varargin)
%function iundist = undistort_image(im, cam, varargin)
%
% im is the original distorted image file (grayscale)
% cam is the camera structure with calibration parameters
%
% varargin{1} is the rectangular image section [x y w h] to be undistorted
% the rectangle is formed from point x,y with width (w) and height (h).
% if not specified the whole image is undistorted
%
% iundist is the undistorted image (or section)
%
% ref: http://www.vision.caltech.edu/bouguetj/calib_doc

% convert image to double
im=double(im);

[h w]=size(im);
if(nargin>2)
    % if a section of image is specified then use that
    rect=varargin{1};
    x=rect(1); y=rect(2);
    wr=rect(3); hr=rect(4);
else
    x=0; y=0;
    wr=w; hr=h;
end

% initialize the undistorted image
iundist=ones(hr,wr);    

% get pixels in the image (or section)
% this is to find where is the portion in real 3D space that must be
% distorted
[xp yp]=meshgrid(x+1:wr+x, y+1:hr+y);

ikm=inv(cam.km);

% find the corresponding values in non-pixel units (mm for example)
x_arr=ikm(1,1)*xp + ikm(1,2)*yp + ikm(1,3);
y_arr=ikm(2,1)*xp + ikm(2,2)*yp + ikm(2,3);

% Apply distortion to them. This gives us a region to apply distortion to.
% We can apply more coefficients and tangential too
r_arr2=x_arr.^2+y_arr.^2;
xd=x_arr.*(1+cam.kc1*r_arr2+cam.kc2*r_arr2.^2);
yd=y_arr.*(1+cam.kc1*r_arr2+cam.kc2*r_arr2.^2);

% back to pixels again
xp_d=xd*cam.km(1,1) + yd*cam.km(1,2) + cam.km(1,3);
yp_d=xd*cam.km(2,1) + yd*cam.km(2,2) + cam.km(2,3);

% now for every distorted [xp_d yp_d], the undistorted one is [xp yp]
px_df=floor(xp_d(:));
py_df=floor(yp_d(:));

px_dc=ceil(xp_d(:));
py_dc=ceil(yp_d(:));

% find points that lie within the limits (after distortion points may go
% out of the region or image)
inds=find((px_df<=wr+x) & (px_df>x) & (py_df<=hr+y ) & (py_df>y));

% update the values so that only those are stored that are within the
% region.
px_df=px_df(inds);
py_df=py_df(inds);
px_dc=px_dc(inds);
py_dc=py_dc(inds);

% find the offset from closest four pixels
ldiff_xp=xp_d(inds)-px_df;
ldiff_yp=yp_d(inds)-py_df;
udiff_xp=1-ldiff_xp;
udiff_yp=1-ldiff_yp;

% pick values from the actual image
pval_lu=im((px_df-1)*h+py_df);
pval_ru=im((px_dc-1)*h+py_df);
pval_ld=im((px_df-1)*h+py_dc);
pval_rd=im((px_dc-1)*h+py_dc);

% find the indices of points that are valid. these are simply the points in
% a 2D matrix that has been reshaped into a column vector. The u,v position
% in 2D is simply (v-v0)*h + u
updated_inds=(xp(inds)-x-1)*hr+yp(inds)-y;

iundist(updated_inds)=(pval_lu).*udiff_xp.*udiff_yp + (pval_ru).*ldiff_xp.*udiff_yp + ...
        (pval_ld).*udiff_xp.*ldiff_yp + (pval_rd).*ldiff_xp.*ldiff_yp;

    
