function wts = p_mq_velocity(Z, tp, cam, Xi, optTrack)
%function wts = p_mq_velocity(Z, tp, cam, Xi, optTrack)

sigma_ep=optTrack.auto.sigma_ep;
te=optTrack.te;

r=tp(Xi.ri,:);
rdot=tp(Xi.rdi,:);%.*(ones(3,1)*tp(Xi.s,:));
% rddot=tp(Xi.rddi,:); 

% r1=r-rdot*te/2-rddot*te^2/4;
% r2=r+rdot*te/2+rddot*te^2/4;

r1=r-rdot*te/2;
r2=r+rdot*te/2;

e1=w2cam(r1,cam);
e2=w2cam(r2,cam);


wts=normal2(e1,Z.ep(:,1), sigma_ep)'.*normal2(e2,Z.ep(:,2), sigma_ep)' + ...
    normal2(e1,Z.ep(:,2), sigma_ep)'.*normal2(e2,Z.ep(:,1), sigma_ep)';


% wts=normpdf(e1(1,:), Z.h(1), sigma_ep(1)).* normpdf(e1(2,:), Z.h(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.t(1), sigma_ep(1)).* normpdf(e2(2,:), Z.t(2), sigma_ep(2));

% wts=normpdf(e1(1,:), Z.h(1), sigma_ep(1)).* normpdf(e1(2,:), Z.h(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.t(1), sigma_ep(1)).* normpdf(e2(2,:), Z.t(2), sigma_ep(2)) + ...
%     normpdf(e1(1,:), Z.t(1), sigma_ep(1)).* normpdf(e1(2,:), Z.t(2), sigma_ep(2)).* ...
%     normpdf(e2(1,:), Z.h(1), sigma_ep(1)).* normpdf(e2(2,:), Z.h(2), sigma_ep(2));
