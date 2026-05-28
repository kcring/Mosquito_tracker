function [Xh P]=initialize1(k, zt, t, Zk, Xh, P, Xh_, P_,  MHT, Kalman, Xi)
global gf;
global idl;

k1=min(k, MHT.nscanback+1); % to search the temporary state space

t=id2tbr(t, gf);
[r c]=getind(Xi.nX, k1, t, 1:Xi.nX, Xi.nX);

% if its a new target it should not be assigned to multiple
% measurements
if numel(zt) > 1
    error('a target cannot be assigned to more than one measurement');
end
[Xh(r,k1,1), P(r,c, 1)]= initTarget(Zk(zt,:), Xi);