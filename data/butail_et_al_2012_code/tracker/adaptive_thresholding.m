function imstream = adaptive_thresholding(p, cam, imstream, cc, Xi, optTrack)
%function imstream = adaptive_thresholding(t_ind, Xh_, cam, imstream, cc, Xi, optTrack)

global gdebug

bt=imstream.binary_t;
foundit=0;
bbox=[1 1 size(imstream.imgarr,2), size(imstream.imgarr,1)];

% ------- stage 1 create an image with lower bt and list the measurements
if isempty(imstream.Zk2)
    imstream.Zk2=getZ(imstream, cc, optTrack, 2*bt/3, bbox, imstream.br);
end

% ---------- stage 2 search for the closest measurement and if it's still
% not there then don't do anything (to avoid double counting measurements)
[val, idx]=min(geticov(p, Xi, imstream.Zk2, cam, optTrack));
% don't change gate size here because it affects in other parts of the code
% extra measurements will be included
if val < optTrack.auto.gatesize
    imstream.Zk=[imstream.Zk, imstream.Zk2(idx)];
    foundit=1;
    if gdebug
        Zplot=cat(2,imstream.Zk2(idx).pixel_list);
        figure(cc+2); gcf; hold on;
        plot(Zplot(:,2), Zplot(:,1), 'c.', 'markersize', 2);
    end
end


if foundit
      fprintf('*)..');
else
    fprintf(' )..');
end
