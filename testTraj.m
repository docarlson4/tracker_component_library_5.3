clc
clear
close all

cd(fileparts(which(mfilename)))

restoredefaultpath

addpath C:/Users/dougo/OneDrive/Documents/MATLAB/
% Constants

addpath(genpath("C:/Users/dougo/source/repos/tracking/tracker_component_library_5.3/"))
Constants;

dbstop if error

% Doug Carlson (doug.o.carlson@gmail.com), 2026-05-24 17:24

%%
traj_obj = Tracker.Trajectory()

s = traj_obj.s
pos = s(1:3,:);

figure
plot(pos(1,:), pos(2,:), '.-')
axis equal
