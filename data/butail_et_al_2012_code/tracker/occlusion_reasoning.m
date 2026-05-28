function [occ n_occ] = occlusion_reasoning(as)
%function [oc1 n_occ] = occlusion_reasoning(as)
% as is the nc x nt matrix where nc is the number of cameras and nt is the
% number of targets tracked
% each element of as is an assignment id of the measurement. This function
% will assign repeated entries in as to occlusions. 
%
% n_occ is the number of occlusions
% oc1 is a structure with size 1 x n_occ
%
%

% initialize
nt=size(as,2);
nc=size(as,1);
k=1;
O1=zeros(nt);
occ=[];

for c = 1:nc % for each camera
	O=ones(nt);
	for i = 1:nt
		for j = (i+1):nt
			% a target can only be assigned to one occlusion
			if isempty(find(O(:,j)==0)) && as(c,i)
				O(i,j)=as(c,i)-as(c,j);
			end
		end
	end

	O(O~=0)=1;

	% now count the occlusions overall
	for i = 1:nt
		if sum(O(i,:))~=nt
			oc1(k).t=[i, find(O(i,:)==0)];
            % O1 has a 1 for each target in an occlusion. The size of O1 is
            % number of occlusions x number of targets
			O1(k,oc1(k).t)=1; 
			k=k+1;
		end
	end
end

% find common occlusions across cameras and union
sr=sum(O1);
for i=1:nt
    if sr(i)>1
        or=find(O1(:,i)==1); % find rows that have the same target in multiple occlusions
        O1(or(1),:)=sum(O1(or,:));
        O1(or(2:end),:)=0; % remove the rest of the occlusions to avoid double counting
        O1(O1>1)=1;
        sr=sum(O1);
    end
end

% now count the occlusions overall
k=1;
for i = 1:nt
    nzi=find(O1(i,:)~=0);
    if ~isempty(nzi)
        occ(k).t=nzi;
        k=k+1;
    end
end

% noc
n_occ=size(occ,2);


% print
if n_occ
    fprintf('Occlusion detected. Occluded targets are:\n');
    for jj=1:n_occ
        fprintf('%d...',occ(jj).t);
        fprintf('\n');
    end
end