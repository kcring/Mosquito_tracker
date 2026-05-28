function Xi=strXi


% state space 
Xi.ri=1:3; % head position
Xi.rdi=4:6; % heading


% nX=0;
% flds=fieldnames(Xi);
% for jj=1:length(flds)
%     nX=nX+numel(Xi.(flds{jj}));
% end
Xi.nX=6;

Xi.cX=[Xi.nX,1];
