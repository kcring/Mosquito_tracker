function [r_sub c_sub] = getind(nx, k, mqid, sv_ind, ncol)
%function [r_sub c_sub] = getind(nx, k, mqid, sv, ncol)
%
% nx: is the number of elements in the state vector
% k: time-step
% mqid: target id
% sv_ind: rows that belong to the value in state vector. for e.g. 1:3 for
% position, 4:6 for velocity
% ncol: number of columns to get
%
% r_sub, c_sub are the row and column index. get the actual values by
% running Xh(r_sub, c_sub)

r_sub=nx*(mqid-1)+sv_ind(1):nx*(mqid-1)+sv_ind(end);
c_sub=ncol*(k-1)+1:ncol*k;