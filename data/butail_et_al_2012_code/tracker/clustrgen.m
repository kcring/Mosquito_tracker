function [C vm2]=clustrgen(vm)
%function C=clustrgen(vm)
% 
% vm as generated from genvmat
% with first column all ones

[nz, nt]=size(vm);
nt=nt-1;

C=strClstr(nz);
for ii=1:nz
    C(ii).z_id(1)=ii;
    t_id=find(vm(ii,2:end));
    C(ii).t_id(1:numel(t_id))=t_id;
end

% combine clusters with common measurements and targets
for ii=1:nz
    for jj=1:nz
        if ii ~=jj
            common=intersect(C(ii).t_id, C(jj).t_id);
            common=common(common~=0);
            if ~isempty(common)
                C(ii).t_id=union(C(ii).t_id, C(jj).t_id);
                C(ii).z_id=union(C(ii).z_id, C(jj).z_id);
                C(jj).t_id=[];
                C(jj).z_id=[];
            end
        end
    end
end

vm2=vm(:,2:end);
for cl=1:size(C,2)
    t=C(cl).t_id;
    z=C(cl).z_id;
    vm2(z,t)=(vm2(z,t)~=0)*cl;
    C(cl).vm=[ones(numel(z),1), vm(z,t)];
end



    