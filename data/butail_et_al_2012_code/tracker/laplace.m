function X=laplace(r,c, mu, b)
% function X=laplace(r,c, mu, b)
% sample from a laplace distribution

U=-.5+rand(r,c);

if nargin ==2
    mu=0; b=1;
end
   

X=mu-b*sign(U).*log(1-abs(U));
