function scan = normalize_prob(k, scan)
        
% ---- Normalize
psum=sum(cat(1,scan(k).hyp.prob));
for hh=1:size(scan(k).hyp,2)
    scan(k).hyp(hh).prob=scan(k).hyp(hh).prob/psum;
end