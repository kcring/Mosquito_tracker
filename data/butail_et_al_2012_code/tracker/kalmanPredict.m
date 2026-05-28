function [Xh_, P_] = kalmanPredict(Xh, P, Klmn)
%function [Xh_, P_] = kalmanPredict(Xh, P, Klmn)
%
% Klmn.F is nxn transition matrix
% Klmn.Q is nxn disturbance matrix


% state
Xh_= Klmn.F*Xh;

% covariance
P_ = Klmn.F*P*Klmn.F' + Klmn.Q;