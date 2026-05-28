function [handles, ind] = t1_triangulate(handles)
% localize the 3D position using the two clicks
jj=handles.mqid;

m(:,1)=handles.frame(handles.k).mq(jj).c(:,1);
m(:,2)=handles.frame(handles.k).mq(jj).c(:,2);


if(handles.backtrack)
    kt=handles.k+1;
else
    kt=handles.k-1;
end

switch handles.optTrack.trackone.alg
    case 'mpf'
        [r, c]=getind(handles.Xi.nX, kt, jj, 1:handles.Xi.nX, 1);
        if(size(handles.Xh,1)>=r(handles.Xi.ri(1)))

            if(handles.Xh(r(1), kt))
                rprev=handles.Xh(r(handles.Xi.ri),kt);
            else
                rprev=[100,100,900]';
            end
        else
            rprev=[100,100,900]';
        end
        % set optimization options
%         options = optimset('Display','off','TolFun',1e-6);
        
%         [fval, lsr]=localize_lsq(rprev, m, handles.get_cam_calib, [1,2]',...
%                                     [-1000, -1000, -50]',[1500, 1500, 4200]');

        [lsr, fval]=lsTriangulate(m, handles.cams);
end

[r, c]=getind(handles.Xi.nX, handles.k, jj, 1:handles.Xi.nX, 1);
handles.Xh(r(handles.Xi.ri), c)=lsr';

if (handles.backtrack)
    if(handles.Xh(r(1), c+1))
        handles.Xh(r(handles.Xi.rdi),c)=(handles.Xh(r(handles.Xi.ri),c+1)-handles.Xh(r(handles.Xi.ri),c))/handles.optTrack.dt;
     end
else
    if(handles.Xh(r(1), c-1))
        handles.Xh(r(handles.Xi.rdi),c)=(handles.Xh(r(handles.Xi.ri),c)-handles.Xh(r(handles.Xi.ri),c-1))/handles.optTrack.dt;
    end
end


% handles.Xh(r(handles.Xi.fi), c)=[-1; 2];

% call splice tracks after every update
% change 9/5/2018 Sachit - add override to the condition
if isfield(handles, 'Xh0') && ~get(handles.cb_override_existing_track, 'Value')
    [handles.Xh, ind]= t1_splice_tracks(handles.Xh, handles.Xh0, ...
            handles.mqid, handles.k, handles.backtrack, handles.Xi, handles.Xi0, handles.optTrack.trackone.threed_t);
else
    ind=0;
end



