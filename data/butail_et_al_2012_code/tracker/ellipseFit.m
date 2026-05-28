function e=ellipseFit(pts)

x=pts(1,:);
y=pts(2,:);

e.center=mean([x;y], 2);
[V D]=eig(cov([x',y']));

if D(2,2)>D(1,1)
    majora=V(:,2);
    length=D(2,2);
else
    majora=V(:,1);
    length=D(1,1);
end


e.majora=majora*sqrt(length)*3;