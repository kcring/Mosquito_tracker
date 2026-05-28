function [S, cost]=murtykbest(costmat, munkres, k)
% function S=murtykbest(costmat, munkres, k)

% reference:
% [1] I. J. Cox and S. L. Hingorani, “An efficient implementation of Reid’s 
% multiple hypothesis tracking algorithm and its evaluation for the purpose 
% of visual tracking,” 
% IEEE Trans. on Pattern Analysis and Machine Intelligence, 
% vol. 18, no. 2, pp. 138-150, 1996.

% number of measurements since we use this in MHT
[nz, nt]=size(costmat);
S=zeros(k, nz);
cost=zeros(k,1);

P=strP(nz,nt);
P(1).costmat=costmat;

% the first solution is from hungarian method
[a, P(1).cost]=munkres(costmat');

% this function gets the assignment
P(1).S=getS(a);

for ii=1:k
    % get the problem with minimum cost
    [val, idx]=min(cat(1,P.cost));
    S(ii,:)=P(idx).S;
    if ii>1 && ~sum(abs(S(ii-1,:)-S(ii,:)))
        S=S(1:ii-1,:);
        break
    end
       
    
    cost(ii)=P(idx).cost;
    P1=P(idx);
    % remove P(idx) from the list of P/S pairs by raising the cost very
    % high so that it is ignored in the next iteration
    P(idx).cost=10^100;
    for jj=1:nz
        P_=P1;
        % remove the assignment from P_ by raising the costmat to a high
        % value 4.4.2
        if P_.S(jj)
            P_.costmat(jj,P_.S(jj))=10^4;
            [a, P_.cost]=munkres(P_.costmat');
            P_.S=getS(a);
            % add P_ to existing problems
            P(end+1)=P_;

            % remove all t and z associations
            P1.costmat(jj, 1:end~=P_.S(jj))=10^4;
            if P_.S(jj)
                P1.costmat(1:end~=jj, P_.S(jj))=10^4;
            end
        end
    end
end
ss=size(S,1);
ntid=max(S(:));
for jj=1:ss
    zi=find(S(jj,:)==0);
    for zz=1:numel(zi)
        ss=ss+1;
        S(ss,:)=S(jj,:);
        S(ss,zi(zz))=ntid+1;
        ntid=ntid+1;
    end
end

        

function P=strP(nz,nt)

P=struct('costmat', zeros(nz,nt), ...
         'cost', 0, ...
         'S', zeros(1,nz));

function assign = getS(a)

[val idx]=max(a,[],1);
assign=idx.*sum(a,1);