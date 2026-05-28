function vec = normv(vec)
%function vec = normv(vec)
%
% to normalize a matrix of vectors
% vec is [2 x n] or [3 x n]

dim=size(vec,1);

vec = vec./(ones(dim,1)*sqrt(sum(vec.^2)));