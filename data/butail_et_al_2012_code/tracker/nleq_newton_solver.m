% Newton's method for nonlinear equations 
% Uses linear search...
%
% fx is equation^2 
% Jx is Jacobian
%
% example: to solve x^2-1==0 run
% nleq_newton_solver(.75, @(x) (x^2-1)^2, @(x) (2*(x^2-1)*2*x), .0001)   
%
% vector example: (to solve two equations x^2-1==0, x^3-3==0, run
% nleq_newton_solver([.75 1]',  @(x) [(x(1)^2-1)^2; (x(2)^3-3)^2], @(x)
% [(2*(x(1)^2-1)*2*x(1)), 0; 0, (2*(x(2)^3-3)*3*x(2))], .0001)
%
% S. Butail 2/10/2011

function sol = nleq_newton_solver(x0, fx, Jx, threshold)

xiter=x0;
iter=1;
while norm(fx(xiter)) > threshold && iter < 50
    step=-Jx(xiter)\fx(xiter);
    xiter=xiter+step;
    iter=iter+1;
%     fprintf('fx(xiter)=%.3f\n', fx(xiter));
end

sol=xiter;