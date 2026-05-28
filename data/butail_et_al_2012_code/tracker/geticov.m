function icov=geticov(p, Xi, Zk, cam, optTrack)

pzh_=w2cam(p(Xi.ri,:), cam);
S=cov(pzh_');
zh_=postest(pzh_,ones(optTrack.auto.Np,1),1);
icov=zeros(size(Zk,2),1);
for jj=1:size(Zk,2)
    icov(jj)=(zh_-Zk(jj).u)'/S*(zh_-Zk(jj).u);
end
if det(S) < 10^-5
    %keyboard;
    fprintf('[!] Badly scaled innovation matrix... \n');
    % catching badly scaled matrix
end    


% % alternate testing
% for jj=1:size(Zk,2)
%     wts=p_lfn(Zk(jj), p, cam, Xi, optTrack);
%     pz(jj)=sum(wts)/numel(wts);
% end
% 
% pz=pz/max(pz);
% fprintf('--icov----pz---\n');
% for jj=1:size(Zk,2)
% 	fprintf('%.2f   %.2f\n', icov(jj), pz(jj));
% end
