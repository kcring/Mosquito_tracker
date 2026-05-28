function scan = hypotheses_reduction(k, scan, MHT)

switch MHT.reduction_strategy
    case 'bestsort'
        [val idx]=sort(cat(1,scan(k).hyp.prob), 'descend');
        if numel(val) >= MHT.sortbest
            scan(k).hyp=scan(k).hyp(idx(1:MHT.sortbest));
            scan=normalize_prob(k, scan);
        end
        
    otherwise
        
end