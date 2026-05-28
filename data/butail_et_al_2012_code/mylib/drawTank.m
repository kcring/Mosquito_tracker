function drawTank(th, tw)
%function drawTank(th, tw)
%
% th, tw are tank height and width

elm=50;
[Cx Cy Cz]=cylinder(tw/2, elm);
    
gcf;
surf(Cx,Cy,Cz*th, 'FaceColor', 'Blue', 'EdgeColor', 'none', 'FaceAlpha', 0.08);hold on;

% base
theta=linspace(0,2*pi, elm);
xb=tw/2*cos(theta);
yb=tw/2*sin(theta);
[Xbase Ybase]=meshgrid(xb, yb);
patch(Xbase', Ybase, ones(elm)*0, 'FaceColor', 'Blue', 'EdgeColor', 'none', 'FaceAlpha', 0.005);