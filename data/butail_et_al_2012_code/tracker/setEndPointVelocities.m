function Zk = setEndPointVelocities(Zk)

global gdebug

x=Zk.pixel_list(:,2);
y=Zk.pixel_list(:,1);
rW=[x';y'];
xaxis=normv(Zk.v);
yaxis=[-xaxis(2); xaxis(1)];
wTb=[xaxis, yaxis, Zk.u;
        0   0       1];
bTw=inv(wTb);    
rB=tra2b(rW,bTw);
xcB=rB(1,:);

p=lineFit(xcB',rB(2,:)');

ycB=p(1)+p(2)*xcB;

[val idx]=sort(xcB);
xcB=xcB(idx); ycB=ycB(idx);

[val e1]=min(xcB); [val e2]=max(xcB);
rcW=tra2b([xcB; ycB], wTb);
for dd=1:2
    rcW(dd,:)=sma(rcW(dd,:), 5);
end
rcv=normv(mydiff(rcW,1));
Zk.length=sum(sqrt(sum(rcv.^2)));
% Zk.speed=Zk.length/(1/40); % pixels/s

Zk.ep(:,1)=rcW(:,e1);
Zk.ep(:,2)=rcW(:,e2);

% update the midpoint
Zk.u=tra2b([0, p(1)]', wTb);

% Zk.ep=rcW(:,[e1, e2]);
% Zk.epv=rcv(:,[e1, e2]);
% Zk.epv(:,1)=-Zk.epv(:,1);
% Zk.vel=Zk.epv*Zk.speed; % pixels/s
% Zk.accel=[Zk.vel(:,1)-Zk.vel(:,2), Zk.vel(:,2)-Zk.vel(:,1)]; % pixel/s^2

if  gdebug % debug
    gca;
%     pt=(-12:12); 
%     xline=wTb(1:2,3)*ones(1,numel(pt))+wTb(1:2,1)*pt;
%     pt=(-4:4);
%     yline=wTb(1:2,3)*ones(1,numel(pt))+wTb(1:2,2)*pt;
%     plot(xline(1,:), xline(2,:), 'Color', ones(1,3)*.5, 'linewidth', 1);
%     plot(yline(1,:), yline(2,:), 'Color', ones(1,3)*.5, 'linewidth', 1);
%     plot(rcW(1,:), rcW(2,:), 'Color', ones(1,3), 'LineWidth', 1.5);
    plot(Zk.u(1), Zk.u(2), 'ro', 'markersize', 6);
    plot(Zk.ep(1,1), Zk.ep(2,1), 'go', 'markersize', 4);
    plot(Zk.ep(1,2), Zk.ep(2,2), 'go', 'markersize', 4);
end

if 0
    figure(10); gcf; 
    
    plot(xcB, ycB-p(1), 'Color', [.5 .5 .5]);
    hold on;
end