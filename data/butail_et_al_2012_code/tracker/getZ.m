function [Z, fg]=getZ(imstream, camid, optTrack, varargin)
%function [Z fg]=getZ(imstream, camid, optTrack, varargin)
% gets measurements...

global gdebug

[h, w]=size(imstream.imgarr(:,:,1));

if nargin>3
    imstream.binary_t=varargin{1};
    
    % --------- we can select a specific region instead
    bbox=varargin{2};
    umin=max(bbox(2), imstream.roi(2)); % rows
    umax=min(bbox(2)+bbox(4), imstream.roi(4));
    vmin=max(bbox(1), imstream.roi(1)); % columns
    vmax=min(bbox(1)+bbox(3), imstream.roi(3));
    
    % for a slow mosquito br is large
    % for a fast mosquito br is small
    imstream.br=varargin{3};
    
    imstream.imgarr=imstream.imgarr(umin:umax,vmin:vmax,:);
end

if(optTrack.fg_is_dark)
    bg=max(imstream.imgarr(:,:,optTrack.br0+1-imstream.br:optTrack.br0+1+imstream.br), [], 3);
else
    bg=min(imstream.imgarr(:,:,optTrack.br0+1-imstream.br:optTrack.br0+1+imstream.br), [], 3);
end

% start hack
idx=ones(1,100);
while numel(idx) > 20 && camid > 0

    fg=imsubtract(bg, imstream.imgarr(:,:,optTrack.br0+1));
    fg=fg/imstream.bitval;


    fg=im2bw(fg, imstream.binary_t);

    if nargin ==3
        fg=setRoi(fg, imstream.roi);
    end
    
    if strcmp(version('-release'), '2007a')
        Li=bwlabel(fg);
    else
        Li=logical(fg);
    end

    
    si = regionprops(Li,'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                        'Orientation', 'PixelIdxList');

%     if size(si,1) > 100
%         error('more than 100 targets detected, check your thresholds!');
%     end                 

    %--------------- speeding it up (calling regionprops only once)...
    area = [si.Area];
    idx = find(area >= imstream.area_t(1) & area < imstream.area_t(2));
    imstream.binary_t=imstream.binary_t*1.05;
end
% this allows the initialization to go as before
if camid < 0
    fg=imsubtract(bg, imstream.imgarr(:,:,optTrack.br0+1));
    fg=fg/imstream.bitval;


    fg=im2bw(fg, imstream.binary_t);

    if nargin ==3
        fg=setRoi(fg, imstream.roi);
    end
    
    if strcmp(version('-release'), '2007a')
        Li=bwlabel(fg);
    else
        Li=logical(fg);
    end


    si = regionprops(Li,'Centroid', 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
                        'Orientation', 'PixelIdxList');
               

    %--------------- speeding it up (calling regionprops only once)...
    area = [si.Area];
    idx = find(area >= imstream.area_t(1) & area < imstream.area_t(2));
end
% end hack

new_fg=zeros*fg;
new_fg(cat(1,si(idx).PixelIdxList))=1;

centroids=cat(1,si(idx).Centroid);
if nargin>3 && ~isempty(centroids)
    centroids(:,1)=centroids(:,1)+vmin;
    centroids(:,2)=centroids(:,2)+umin;
end

majora=cat(1,si(idx).MajorAxisLength);
minora=cat(1,si(idx).MinorAxisLength);
orientation=cat(1,si(idx).Orientation);
area=cat(1,si(idx).Area);

[hi, wi]=size(imstream.imgarr(:,:,1));

if nargin>3
    tmp_fg=zeros(h,w);
    tmp_fg(umin:umax,vmin:vmax)=new_fg;
    fg=tmp_fg;
else
    fg=new_fg;
end


if  gdebug % for debug 
    figure(10); gcf;
    j=(uint16(imstream.imgarr(:,:,optTrack.br0+1))); %histeq
%     imshow(j); hold on;
    imshow((1-logical(fg))*255); hold on;
    
end

Z=strZ(numel(idx));
for ii=1:numel(idx)
    Z(ii).u=centroids(ii,:)';
    % ------- orientation vector
%     jj=[majora(ii)*cosd(orientation(ii)); -majora(ii)*sind(orientation(ii))]/2;
%     Z(ii).ep=[Z(ii).u+jj Z(ii).u-jj];
    Z(ii).v=[majora(ii).*cosd(orientation(ii)), -majora(ii).*sind(orientation(ii))]';
    % ------- frame transformation bTw = [cost +sint; -sint cost],
    % therefore wTb = bTw'
    Z(ii).sigma=sqrt(abs([majora(ii).*cosd(orientation(ii)) - minora(ii).*sind(orientation(ii)); 
                          majora(ii).*sind(orientation(ii)) + minora(ii).*cosd(orientation(ii))]));
    Z(ii).area=area(ii);
    [i1, i2]=ind2sub([hi, wi], si(idx(ii)).PixelIdxList);

    if nargin>3
        Z(ii).pixel_list=[i1+umin-1, i2+vmin-1];
    else
        Z(ii).pixel_list=[i1,i2];
    end
    
    Z(ii)=setEndPointVelocities(Z(ii));
    % using the index of position as the id
    Z(ii).id=sub2ind([h,w], ceil(Z(ii).u(2)), ceil(Z(ii).u(1)));
end



