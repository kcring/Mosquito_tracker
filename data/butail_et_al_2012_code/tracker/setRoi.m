function fg = setRoi(fg, roi)
%function fg = setRoi(fg, roi)

% remove everything except roi
fg(1:roi(2),:)=0;
fg(roi(2)+roi(4):end,:)=0;
fg(:,1:roi(1))=0;
fg(:,roi(1)+roi(3):end)=0;