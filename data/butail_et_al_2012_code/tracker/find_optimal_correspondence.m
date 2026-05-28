function [z1h, z2h]=find_optimal_correspondence(z1,z2,F)

x=z1(1); y=z1(2);
xp=z2(1); yp=z2(2);

% (i)
T=[1 0 -x; 0 1 -y; 0 0 1];
Tp=[1 0 -xp; 0 1 -yp; 0 0 1];


% (ii)
F=Tp'\F/(T);


% (iii)
e=null(F); f=e(3);
e = e(1:2)/norm(e(1:2));
ep=null(F'); fp=ep(3);

ep=ep(1:2)/norm(ep(1:2));


% (iv)
R=[e(1) e(2) 0; -e(2) e(1) 0; 0 0 1];
Rp=[ep(1) ep(2) 0; -ep(2) ep(1) 0; 0 0 1];

% (v)
F=Rp*F*R';

% (vi)
a=F(2,2); b=F(2,3); c=F(3,2); d=F(3,3);

% (vii)
% gt=t*((a*t+b)^2+fp^2*(c*t+d)^2)^2-(a*d-b*c)*(1+f^2*t^2)^2*(a*t+b)*(c*t+d)
% ;

gt(7)=b*d*(b*c - a*d);
gt(6)=b^4 + b*c*(b*c - a*d) + a*d*(b*c - a*d) + 2*b^2*d^2*fp^2 + d^4*fp^4;
gt(5)=4*a*b^3 + a*c*(b*c - a*d) - 2*b*d*(b*c - a*d)*f^2 + 4*b^2*c*d*fp^2 + ... 
    4*a*b*d^2*fp^2 + 4*c*d^3*fp^4;
gt(4)=6*a^2*b^2 - 2*b*c*(b*c - a*d)*f^2 - 2*a*d*(b*c - a*d)*f^2 + ...
    2*b^2*c^2*fp^2 + 8*a*b*c*d*fp^2 + 2*a^2*d^2*fp^2 + 6*c^2*d^2*fp^4;
gt(3)=4*a^3*b - 2*a*c*(b*c - a*d)*f^2 + b*d*(b*c - a*d)*f^4 + ...
    4*a*b*c^2*fp^2 + 4*a^2*c*d*fp^2 + 4*c^3*d*fp^4;
gt(2)= a^4 + b*c*(b*c - a*d)*f^4 + a*d*(b*c - a*d)*f^4 + 2*a^2*c^2*fp^2 + ...
    c^4*fp^4;
gt(1)=a*c*(b*c - a*d)*f^4; 

rgt=roots(gt);
t=real(rgt);

% (viii)
st=t.^2./(1+f^2*t.^2)+(c*t+d).^2./((a*t+b).^2+fp^2*(c*t+d).^2);

t(7)=inf;
st(7)=1/f^2+c^2/(a^2+fp^2*c^2);  % at t=inf;

[val idx]=min(st);

tmin=t(idx);

if tmin > 10000
    fprintf('[!] tmin at inf..\n');
end
t=tmin;

% (ix)
l=[t*f; 1; -t];
lp=[-fp*(c*t+d); a*t+b; c*t+d];

z1h=closest_to_origin(l);
z2h=closest_to_origin(lp);

% (x)
z1h=T\R'*z1h;
z2h=Tp\Rp'*z2h;

% extra step to make homogeneous
z1h=z1h/z1h(3);
z2h=z2h/z2h(3);


function zh=closest_to_origin(l)

lamda=l(1);
mu=l(2);
nu=l(3);

zh=[-lamda*nu; -mu*nu; lamda^2+mu^2];
