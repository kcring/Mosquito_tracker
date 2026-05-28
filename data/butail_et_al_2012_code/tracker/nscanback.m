function scan = nscanback(k, scan, MHT)

% scan back is performed in any case
if k-MHT.nscanback > 1 && ~isempty(scan(k-MHT.nscanback).hyp);
    [val idx]=max(cat(1,scan(k).hyp.prob));
    parent=scan(k).hyp(idx).parent;
    % only need this for nscanback > 1
    jj=k-1;
    while jj>k-MHT.nscanback
        parent=scan(jj).hyp(parent).parent;
        jj=jj-1;
    end
    
    % start pruning hypothesis from k-MHT.scanback onwards 
    % keeping only the one that has the parent from
    try
    % first the only hypothesis surviving at nscanback is the 
    % great grandparent only    
    scan(k-MHT.nscanback).hyp=scan(k-MHT.nscanback).hyp(parent);

    catch
        keyboard
    end
    for jj=k-MHT.nscanback+1:k
        % search for all the parents
        parents=cat(1,scan(jj).hyp.parent);
        % look for great grandparent
        idx=find(ismember(parents,parent));
        % update the hypotheses
        scan(jj).hyp=scan(jj).hyp(idx);
        % update the parent id to newparent
        for pp=1:numel(idx)
            scan(jj).hyp(pp).parent=find(parent==scan(jj).hyp(pp).parent);
        end
        
        parent=idx;
        scan=normalize_prob(jj,scan);
    end
    
end