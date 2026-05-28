function m_mqid=show_tracked_mq(Xh, Xi, k, disp)

if ~k,  k=1:size(Xh,2); end


m_flags=Xh(Xi.ri(3):Xi.nX:end,k);
m_ind=find(m_flags~=0);
[m_mqid k]=ind2sub(size(m_flags), m_ind);

m_mqid=unique(m_mqid);

if disp
    if(~isempty(m_mqid))
        fprintf('Total trackone tracks = %d ... \n', numel(m_mqid));
        fprintf('\n');    
        fprintf('id \t start \t end\n');
        fprintf('---------------------\n');
        for mq_ii=m_mqid'
            [r c]=getind(Xi.nX, 1, mq_ii, 1, 1);
            nz_ind=find(Xh(r,:)~=0);
            fprintf('%d \t %d \t %d\n', mq_ii, nz_ind(1), nz_ind(end));
        end
        fprintf('=====================\n');    
    end
end
