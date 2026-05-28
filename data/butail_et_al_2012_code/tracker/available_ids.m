function avid=available_ids(tXh, Xi, maxid)
% assuming that r_3 cannot be zero for a valid mosquito
m_flags=tXh(Xi.ri(3):Xi.cX(1):end,:);
m_ind=find(m_flags~=0);
[m_mqid k]=ind2sub(size(m_flags), m_ind);

avid=min(find(~ismember([1:maxid], m_mqid)));
if isempty(avid)
    avid=size(tXh,1)/Xi.cX(1) + 1;
end