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
pos_plat_geo0 = [15*deg; -75.0*deg; hp];
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

%% Set Sensor Location and Attitude on Board A/C

%% Generate Measurements of Targets as Seen by the Sensor
%   - slant range
%   - direction cosine

%% Impose Measurement Noise

%% Reconstruct Target Geolocation Based On Terrain Constraint
%   - use cubature to generate moments: mean and covariance
%       a. geodetic
%       b. Cartesian (tracker measurement and motion model)

%% Plots
figure
hold on, grid on, box on
plot3(pos_plat_ecf(1,:), pos_plat_ecf(2,:), pos_plat_ecf(3,:), LineWidth=1)

% Load the coastline data
load coastlines

% Create a sphere
a = Constants.WGS84SemiMajorAxis;
b = Constants.WGS84SemiMinorAxis;
[X, Y, Z] = ellipsoid(0,0,0,a,a,b);  % Adjust resolution as needed (100 is medium-res)

surf(X, Y, Z, 'FaceColor', 'none', 'EdgeColor', 'none')  % Create transparent sphere

% Plot the coastlines on the sphere
% Convert coastline lat/lon to 3D Cartesian coordinates
radius = 0.99;  % Slightly less than 1 to avoid overlap with sphere
coast_x = a * radius * cosd(coastlat) .* cosd(coastlon);
coast_y = a * radius * cosd(coastlat) .* sind(coastlon);
coast_z = b * radius * sind(coastlat);
coast = [coast_x, coast_y, coast_z];

% Find points on the far side (where z < 0)
visible = coast * mean(state_plat(1:3,:),2) >= 0;  % Set visibility condition (z >= 0 is near side)

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

