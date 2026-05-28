function prob = p_pos(Z,X,cam,sigma)
%function prob = p_pos(Z,X, cam, sigma)
% 
% Likelihood function P(Z|X) for position
% 
% Z is 2 X 1 measurement matrix of silhouette points
% X is nx X N state matrix (nx is the number of states)
% cam is the camera structure
% sigma is the standard deviation matrix
%
% prob is 1 X N matrix


pix=w2cam(X, cam);

% assuming independence along x and y dimensions
prob=normpdf(pix(1,:), Z(1), sigma(1,1)).*...
     normpdf(pix(2,:), Z(2), sigma(2,2));