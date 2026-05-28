function [optTrack, swarm_boundaries, camid_suffixes]=config
% Use this file to configure tracking variables

% Common options
optTrack.figless=1;
optTrack.te=1/40; % exposure time
optTrack.dt=1/60;
optTrack.bbox=[50,50]; %half height and width of bounding box for searching in trackone
optTrack.fg_is_dark=1;
optTrack.br0=3; % sliding window size = 2*d+1

global gdebug
gdebug=0;

swarm_boundaries=[200 900]; % mm

% camid suffixes
camid_suffixes=['L'; 'R'];

