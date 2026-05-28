function hyparr=strhyp(assignments, parent)

hyp=struct('parent', parent, ...
           'ch', [], ...
           'assignment', [], ...
           'prob', 0);
nh=size(assignments,1);
hyparr(1:nh)=hyp;

for hh=1:nh
    hyparr(hh).assignment=assignments(hh,:);
end
