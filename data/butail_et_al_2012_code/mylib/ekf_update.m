function [Xh P]=ekf_update(Xh_, P_, R, Z, Hfun, Hlinfun, varargin)
%function [Xh P]=ekf_update(Xh_, P_, R, Z, Hfun, Hlinfun, varargin)
%
% extended Kalman filter (single iteration)
%
% NOTE: all vectors are row vectors unless specified
%
% Xh_ is the state estimate at the current time step 
% P_ is the state covariance matrix at the current time step
% R is the measurment noise covariance matrix
% Z is the measurement
% Hfun is the nonlinear measurement function handle
% Hlinfun is the linearized measurement function handle
%
% varargin is for later for JPDA
%
% example call:
% 
% [Xh P]=ekf_update(Xh_, P_, R, Z, @(Xh_) Hfun(Xh_, a, b), @(Xh_) Hlinfun(Xh_, a, b))
% where a, b are other arguments that may be required to pass the
% measurement functions

% get state vector
nx = size(Xh_,1);

% get Hk evaluated at current state
Hk = Hlinfun(Xh_);

% measurement covariance
S = Hk*P_*Hk' + R;

% projected estimate
Zh = Hfun(Xh_);

% innovation
inov = Z-Zh;

% Kalman gain
W = P_*Hk\(S);

[msg msgid]=lastwarn;
if strcmp(msgid, 'MATLAB:illConditionedMatrix') || strcmp(msgid, 'MATLAB:singularMatrix')
    keyboard;
end 

% state update
Xh = Xh_ + W*inov;

% covariance update
P = (eye(nx)-W*Hk)*P_;

