function [Xh_ P_]=discrete_cv(Xh, P, dt, vd)
%function [Xh_ P_]=discrete_cv(Xh, P, dt, vd)
%
% discrete constant velocity motion model
%
% Xh is 6x1 vector with 1:3 position and 4:6 velocity
% P is 6x6 covariance matrix
% dt is the time-step size
% vd is the disturbance coefficient



% motion model (constant velocity)
Fk= [  1,  0,  0, dt,  0,  0
       0,  1,  0,  0, dt,  0
       0,  0,  1,  0,  0, dt
       0,  0,  0,  1,  0,  0
       0,  0,  0,  0,  1,  0
       0,  0,  0,  0,  0,  1];
   
   
Qk= [ 1/3*dt^3*vd^2,             0,             0, 1/2*dt^2*vd^2,             0,             0
                  0, 1/3*dt^3*vd^2,             0,             0, 1/2*dt^2*vd^2,             0
                  0,             0, 1/3*dt^3*vd^2,             0,             0, 1/2*dt^2*vd^2
      1/2*dt^2*vd^2,             0,             0,       dt*vd^2,             0,             0
                  0, 1/2*dt^2*vd^2,             0,             0,       dt*vd^2,             0
                  0,             0, 1/2*dt^2*vd^2,             0,             0,       dt*vd^2];

% update the state estimate
Xh_=Fk*Xh;

% update the covariance matrix for next iteration
P_=Fk*P*Fk' + Qk;              