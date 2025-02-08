clc
clear
close all

cd(fileparts(which(mfilename)))

MyConstants

dbstop if error

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2025-02-07 20:13

% Airborne Radar Simulation for Tracking Ground Targets

%% Outline

% Set time span
Ns = 15;
Ts = 10; % sec
t = (0:Ns-1)*Ts;

% The continuous-time drift function for a 3D, flat-Earth dynamic model is
aDrift = @(x,t) aCircTurn(x);

%% Set A/C kinematics

%   1. initial position (GEO/ECEF)
hp = 15*km;
pos_plat_geo0 = [14.90*deg; -75.0*deg; hp];
pos_plat_ecf0 = ellips2Cart(pos_plat_geo0);

%   2. initial velocity (ENU)
va = 300*mph;
pha = 0*deg;
vel_plat_enu0 = va * [sin(pha); cos(pha); 0];

%   3. propagate platform
% The initial target state with GLOBAL position and LOCAL velocity
state_plat_init = [pos_plat_ecf0; vel_plat_enu0; 0.0];

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
[state_plat] = StateProp3D(state_plat_init, t, aDrift);

pos_plat_ecf = state_plat(1:3,:);
pos_plat_geo = Cart2Ellipse(pos_plat_ecf);
vel_plat_enu = state_plat(4:6,:);

%% Set Target Kinematics

%   1. initial position (GEO/ECEF)
ht = 0*km;
pos_tgt_geo0 = [15*deg; -74.0*deg; ht];
pos_tgt_ecf0 = ellips2Cart(pos_tgt_geo0);

%   2. initial velocity (ENU)
vt = 30*knots;
pha = 45*deg;
vel_tgt_enu0 = vt * [sin(pha); cos(pha); 0];

%   3. propagate target
% The initial target state with GLOBAL position and LOCAL velocity
state_tgt_init = [pos_tgt_ecf0; vel_tgt_enu0; 0.0];

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
[state_tgt] = StateProp3D(state_tgt_init, t, aDrift);

pos_tgt_ecf = state_tgt(1:3,:);
pos_tgt_geo = Cart2Ellipse(pos_tgt_ecf);
vel_tgt_enu = state_tgt(4:6,:);

%% Set Sensor Location and Attitude on Board A/C
plhPoint = pos_plat_geo;
% Port side
az = 90*deg;
% Level flight (A/C not banking)
el = 0*deg;
% Level flight (A/C not pitching)
zRot = 0*deg;
M = zeros(3,3,Ns);
for kS = 1:Ns
    M(:,:,kS) = findRFTransParam(plhPoint(:,kS), az, el, zRot);
end

%% Generate Measurements of Targets as Seen by the Sensor
%   - slant range
%   - direction cosines

zC = pos_tgt_ecf;
useHalfRange = true;
zTx = pos_plat_ecf;
zRx = pos_plat_ecf;
includeW = false;
[z] = Cart2Ruv(zC, useHalfRange, zTx, zRx, M, includeW);
% diag([1/km,1,1])*z

% Slant range (m)
rs = z(1,:);
% Direction cosines in ACS
ua = z(2,:);
va = z(3,:); % generally unknown for ULAs

%% Spectrum
f0 = 3*GHz;
lam = c0/f0;
k0 = 2*pi/lam;

%% Antenna
Ne = 8; % no. of elements
d = lam/2;

%% Impose Measurement Noise
radar_obj = Tracker.RadarReceiver( ...
    "CenterFreq", f0, ...
    "FrameTime", Ts);
sig_rs = radar_obj.RangeUnc;
sig_ua = sqrt(3)/(k0*d*sqrt(Ne*(Ne-1)*(2*Ne-1)) * radar_obj.SNR);
SR = diag([sig_rs, sig_ua]);

%% Reconstruct Target Geolocation Based on Terrain Constraint
%   - use cubature to generate moments: mean and covariance
%       a. geodetic
%       b. Cartesian (tracker measurement and motion model)

%% Plots
% Define the size and origin of the scene
% Latitude of the origin
origin_lat = mean([mean(pos_plat_geo(1,:)), mean(pos_tgt_geo(1,:))]/deg);
% Longitude of the origin
origin_lon = mean([mean(pos_plat_geo(2,:)), mean(pos_tgt_geo(2,:))]/deg) ;

scene_width  = 2; % Width of the scene in degrees
scene_height = 1; % Height of the scene in degrees

% Plot the scene with specified origin and size
figure;

% Plot platform path
geoplot(pos_plat_geo(1,:)/deg, pos_plat_geo(2,:)/deg, 'r-', ...
    'LineWidth', 2);
hold on;

% Plot Target path
geoplot(pos_tgt_geo(1,:)/deg, pos_tgt_geo(2,:)/deg, 'b-', ...
    'LineWidth', 2); 

% Set the limits of the scene
geolimits( ...
    origin_lat + [-1/2,1/2]*scene_height, ...
    origin_lon + [-1/2,1/2]*scene_width);

% Set the base map style
geobasemap satellite

title('Movement of Two Objects in Specified Scene');

% Add a legend
legend('Object 1', 'Object 2');

% Display the plot
hold off;

return
figure
hold on, grid on, box on
plot3(pos_plat_ecf(1,:), pos_plat_ecf(2,:), pos_plat_ecf(3,:), ...
    LineWidth=1, DisplayName="Platform")
plot3(pos_plat_ecf(1,:), pos_plat_ecf(2,:), pos_plat_ecf(3,:), ...
    LineWidth=1, DisplayName="Platform")

% Load the coastline data
load coastlines

% Create a sphere
a = Constants.WGS84SemiMajorAxis;
b = Constants.WGS84SemiMinorAxis;

% Adjust resolution as needed (100 is medium-res)
[X, Y, Z] = ellipsoid(0,0,0,a,a,b);

% Create transparent sphere
surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 'none')

% Plot the coastlines on the sphere
% Convert coastline lat/lon to 3D Cartesian coordinates
radius = 0.99;  % Slightly less than 1 to avoid overlap with sphere
coast_x = a * radius * cosd(coastlat) .* cosd(coastlon);
coast_y = a * radius * cosd(coastlat) .* sind(coastlon);
coast_z = b * radius * sind(coastlat);
coast = [coast_x, coast_y, coast_z];

% Find points on the far side (where z < 0)
%   - set visibility condition (z >= 0 is near side)
visible = coast * mean(state_plat(1:3,:),2) >= 0;

% Break the line by inserting NaNs for far side points
coast_x(~visible) = NaN;
coast_y(~visible) = NaN;
coast_z(~visible) = NaN;

% Plot only visible parts of the coastline
plot3(coast_x, coast_y, coast_z, 'k')

ax = gca;
ax.Clipping = 'off';
view(0,35)
axis equal
axis off

