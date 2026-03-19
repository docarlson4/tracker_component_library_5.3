clc
clear
close all

cd(fileparts(which(mfilename)))

MyConstants

dbstop if error

% Doug Carlson (doug.o.carlson@gmail.com), 2024-11-16 20:29

%% Radar
radar_obj = Tracker.RadarReceiver( ...
    "FrameTime", 10);
Ts = radar_obj.FrameTime;
num_scan = 30;

%% Platform
pos_plat_enu = [0;0;0];
vel_plat_enu = [0;0;0];
stt_plat_enu = [pos_plat_enu;vel_plat_enu];

%% Target
tgt_mot_mdl_obj = Tracker.MotionModel( ...
    "Type", "NCA", ...
    "RevisitTime", Ts, ...
    "SigmaAcc", 0.1);
sDim = tgt_mot_mdl_obj.SpaceDim;
xDim = tgt_mot_mdl_obj.StateDim;

pos_tgt_enu = [-30;20;10]*km;
vel_tgt_enu = [200;0;0];
acc_tgt_enu = [0;0;0];
stt_tgt_enu = [pos_tgt_enu;vel_tgt_enu;acc_tgt_enu] * ones(1,num_scan);

for kS = 2:num_scan
    stt_tgt_enu(:,kS) = tgt_mot_mdl_obj.f(stt_tgt_enu(:,kS-1));
end

%% Measurements
cases = [
    "Range Only"
    "Range-Azimuth-Elevation"
    "RAED"];
chosen_case = 3;

switch cases(chosen_case)
    
    case 'Range Only'
        h = @(s) [
            getRange(s(1:3,:), true, pos_plat_enu, pos_plat_enu);
            getRangeRate(s(1:6,:), true, stt_plat_enu, stt_plat_enu)];
        SR = diag([radar_obj.RangeUnc, radar_obj.RangeRateUnc]);
        zDim = size(SR,1);

    case 'Range-Azimuth-Elevation'
        h = @(s) Cart2Sphere(s(1:3,:), 3, true, pos_plat_enu, pos_plat_enu);
        ElevationUnc = sqrt(pi/4) * 15*deg./sqrt(2*radar_obj.SNR);
        SR = diag([radar_obj.RangeUnc, radar_obj.AzimuthUnc, ElevationUnc]);
        zDim = size(SR,1);

    case 'RAED'
        h = @(s) [[1,0,0;0,1,0] * ...
            Cart2Sphere(s(1:3,:), 3, true, pos_plat_enu, pos_plat_enu);
            20*deg * ones(1,size(s,2));
            getRangeRate(s(1:6,:), true, stt_plat_enu, stt_plat_enu)];
        ElevationUnc = sqrt(pi/4) * 15*deg./sqrt(2*radar_obj.SNR);
        SR = diag([radar_obj.RangeUnc, radar_obj.AzimuthUnc, ...
            ElevationUnc, radar_obj.RangeRateUnc]);
        zDim = size(SR,1);

end

zMeas = h(stt_tgt_enu) + SR * randn(zDim, num_scan);

%% Initialization
x = stt_tgt_enu(:,1) * ones(1,num_scan);
SP = blkdiag(1000*eye(sDim), 10*eye(sDim), 1*eye(sDim));
x(:,1) = x(:,1) + SP*(randn(xDim,1));
SP = repmat(SP, [1,1,num_scan]);

%% Cubature Filter
f = tgt_mot_mdl_obj.f;
SQ = cholSemiDef(tgt_mot_mdl_obj.Q, 'lower');
for kS = 2:num_scan

    % Prediction
    xPrev = x(:,kS-1);
    SPrev = SP(:,:,kS-1);
    [xPred, SPred] = sqrtDiscCubKalPred(xPrev,SPrev,f,SQ);
    
    % Update
    [x(:,kS), SP(:,:,kS)] = sqrtCubKalUpdate(xPred,SPred,zMeas(:,kS),SR,h);

    disp('')
end

%% Plots
figure
hold on, grid on, box on
xlabel("\bf East ")
ylabel("\bf North ")
zlabel("\bf Up ")
axis equal

% php = plot3(pos_plat_enu(1)/km, pos_plat_enu(2)/km, pos_plat_enu(3)/km, ...
%     '^', "LineWidth", 1.5, "DisplayName","Platform");

pht = plot3(stt_tgt_enu(1,:)/km, stt_tgt_enu(2,:)/km, stt_tgt_enu(3,:)/km, ...
    'ko', "LineWidth", 1.5, "DisplayName","Target Truth");

for kS = 1:num_scan
    phe = plot3(x(1,kS)/km, x(2,kS)/km, x(3,kS)/km, ...
        'r.', "LineWidth", 1.5, "DisplayName","Pos Est");
    disp('')
end









