function [Xh P]=nnf_kf(Z, cZ, Xh_, cX, P_, cP, Rk, Hfun, Hlinfun, cost, gate)
%function [Xh P]=nnf_kf(Z, Xh_, P_, Rk, Hk, costfun, gate)
%
% Nearest neighbor filter (Kalman filter update)
%
%
% Input:
%        Z:       m*nz x 1 where m is the size of each observation and n is the
%                   number of observations
%        Xh_:     n*nt x 1 where n is the size of the state vector and nt is the number of targets   
%        P_:      n*nt x n     
%        Rk:      m x m         
%        Hfun:    function, set to Hlinfun if Kalman filter
%        Hlinfun: linearized Hfun, set to @(x) Hk for Kalman filter where
%                       Hk is constant
%        cost:    costfunction based on which distance metric is computed
%        gate:    gate size beyond which measurement is not assigned

nt=size(Xh_,1)/cX(1);
nz=size(Z,1)/cZ(1);
Hk=zeros(cZ(1), cP(1), nt);
Xh=Xh_;
P=P_;

% Create inovation covariance matrices
for t=1:nt
    rt=snipx(cX,t, 1:cP(1));
    Hk(:,:,t) = Hlinfun(Xh_(rt, 1));
    S(:,:,t)=Hk(:,:,t)*P_(snipx(cP,t), :)*Hk(:,:,t)' + Rk;
    zh_(snipx(cZ,t),1)=Hfun(Xh_(rt,1));
end

costmat=zeros(nz,nt);

if gate
    inov=zeros(cZ(1),nt);
    gatevol=zeros(nz,nt); 
end

% evaluate cost
for i=1:nz % measurements
    for j=1:nt % estimated measurements corresponding to targets
        rj=snipx(cZ,j);
        ri=snipx(cZ,i);
        if gate
            inov(:,j)=Z(ri,1)-zh_(rj,1);
            gatevol(i,j)=inov(:,j)'/(S(:,:,j))*inov(:,j);
        end
        costmat(i,j)=cost(zh_(rj,1), Z(ri,1));
    end
end

% as is a horizontal vector with target t assigned measurement as(t)
[cost as]=min(costmat, [], 1);

for t=1:nt
    if gatevol(as(t),t) > gate
        as(t)=0;
    end
end

for t=1:nt

    % compute the gain matrix for each target
    rt=snipx(cX,t,1:cP(1));
    rp=snipx(cP,t);
    W(:,:,t)=P_(rp,:)*Hk(:,:,t)'/S(:,:,t);

    % TBD ...
    if as(t)
        inov=Z(snipx(cZ,as(t)),1)-zh_(snipx(cZ,t),1);
    else
        inov=zeros(cZ);
    end
        
    % update the state estimate using combined innovation
    Xh(rt,:)=Xh_(rt,:)+W(:,:,t)*inov;

    % update state covariance
    P(rp,:)=(eye(cP)-W(:,:,t)*Hk)*P_(rp,:);
end    