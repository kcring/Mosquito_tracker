function t1_save_tracks(handles)

% pathname=fileparts(datfile);

Xh=handles.Xh;
Xi=handles.Xi;
frmlist=handles.frmlist;
cams=handles.cams;
if isfield(handles, 'frame')
    frames=handles.frame;
else
    frames=[];
end

save(handles.datfile1, 'frames', 'Xh', 'Xi', ...
                       'frmlist', 'cams', ...
                      '-append');


% m_mqid=show_tracked_mq(Xh, Xi, 0, 'h', 0);
% 
% if(~isempty(m_mqid))
%     for ii=m_mqid'
%         [r c]=getind(Xi.nX, 1, ii, 1:Xi.nX, 1);
%         fprintf('Saving mq-%.2d.csv ....\n', ii);
%         csvwrite(sprintf('%s/mq-%.2d.csv', pathname, ii), Xh(r,:));
% %         [r c]=getind(Zi.nZ, 1, ii, 1:Zi.nZ, 1);
% %         csvwrite(sprintf('%s/Zmq-%d.csv', pathname, ii), Z(r,:));
%     end
% end
