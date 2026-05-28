function [x y]=draw_paralellogram(a, b, theta, corner, draw_yn)


x=corner(1)+[0, a*cos(theta), a*cos(theta), 0];
y=corner(2)+[0, a*sin(theta), a*sin(theta)+b, b];

if draw_yn
    patch(x,y, [.75 .75 .75], 'FaceAlpha', .1, 'EdgeColor', [.5 .5 .5])
end