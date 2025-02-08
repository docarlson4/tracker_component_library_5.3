clc
clear
close all

cd(fileparts(which(mfilename)))

addpath C:/Users/dougo/OneDrive/Documents/MATLAB/
Constants

dbstop if error

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-11-03 12:33

%% Radar Parameters
radar_obj = Tracker.RadarReceiver( ...
    "Frametime", 6, ...
    "Bandwidth", 2*MHz, ...
    "NumPulse", 10, ...
    "SNRdB", 10);

%% Target Swarm
tgt_obj = Tracker.TargetSwarm( ...
    'RadarObject', radar_obj ...
    )












