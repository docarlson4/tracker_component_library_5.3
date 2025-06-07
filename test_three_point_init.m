clc
clear
close all

cd(fileparts(which(mfilename)))

restoredefaultpath

addpath C:/Users/dougo/OneDrive/Documents/MATLAB/
Constants

addpath(genpath("C:/Users/dougo/source/repos/tracking/tracker_component_library_5.3/"))
Constants;

dbstop if error


% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2025-06-07 09:55

%%
Ts = 12;
maxAcc = 9.8/5; % in m/s/s
maxJrk = 0.1;
tau = 25*Ts;   % maneuver decorrelation time (s)
motion_model = Tracker.MotionModel( ...
    "Type", "NCA", ...
    "SpaceDim", 2, ...
    "RevisitTime", Ts, ...
    "MaxKinVar", maxAcc, ...
    "Tau", tau);
% State motion model method - see displayTracksSwarm.m
motion_model_method = motion_model.Type;

% State Initialization
init_methods = ["One-Point", "Two-Point", "Three-Point"];
state_init = Tracker.StateInitialization( ...
    "Type", "Three-Point", ...
    "VelMin", 100, ...
    "VelMax", 150, ...
    "MotionModel", motion_model);
% State initialization method - see displayTracksSwarm.m
init_method = state_init.Type;


%%
Ns = 4;
tMmt = {[0,0,0];[12,12];[24,24,24];[36,36,36]};
z1 = [0,1,3;2,1,-1]*km;
z2 = [2,3;2,1]*km;
z3 = [2,3,4;4,3,2]*km;
z4 = [4,5,6;4,3,2]*km;
zMmt = {z1;z2;z3;z4};
SRMmt = cell(Ns,1);

figure
hold on, grid on, box on
for kS = 1:Ns
    scatter(zMmt{kS}(1,:),zMmt{kS}(2,:))
    SRMmt{kS} = repmat(eye(2),1,1,numel(tMmt{kS}));
end

for kS = 1:Ns
    % State initialization
    tCur = tMmt{kS};
    zCur = zMmt{kS};
    SRCur = SRMmt{kS};
    [xNew, SNew, lgclNew] = state_init.Initialize(tCur, zCur, SRCur)
    scatter(xNew(1,:), xNew(2,:), '*')
    disp('')
    
end







