function [Xh_ P_]=ekf_predict(Xh, P, Qk, Ffun, Flinfun, varargin)
%function [Xh_ P_]=ekf_predict(Xh, P, Qk, Ffun, Flinfun, varargin)
%
% extended Kalman filter (single iteration)
%
% NOTE: all vectors are row vectors unless specified
%
% Xh is the state estimate at the current time step 
% P is the state covariance matrix at the current time step
% Qk is the disturbance covariance matrix
% Ffun is the nonlinear measurement function handle
% Flinfun is the linearized measurement function handle
%
% varargin is for later for JPDA
%
% example call:
% 
% [Xh_ P_]=ekf_predict(Xh, P, Qk, @(Xh) Ffun(Xh, a, b), @(Xh) Flinfun(Xh, a, b))
% where a, b are other arguments that may be required to pass the
% measurement functions

% state
Xh_= Ffun(Xh);

% covariance
P_ = Flinfun(Xh)*P*Flinfun(Xh)' + Qk;