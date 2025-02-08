clc
clear
close all

cd(fileparts(which(mfilename)))

MyConstants

dbstop if error

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-11-03 14:48

%%

%The latitude and longitude of the Mauna Kea Observatory in Hilo, Hawaii.
phi=35.823*deg;%North latitude in radians.
lambda=-85.470*deg;%East longitude in radians.
latLonStart=[phi;lambda];

%Obtain the initial Cartesian location of the target.
ellipsAlt = 8*km;
rGlobalCart=ellips2Cart([latLonStart;ellipsAlt]);

%The initial velocity in LOCAL coordinates of the target is Mach 1 times
%the direction obtained by solving the indirect geodesic problem. The
%direction is found from trigonometry using the heading (in radians East of
%North).
cSoundSTP=speedOfSoundInAir();
azStartE = 45*deg;
vLocalInit=cSoundSTP*[sin(azStartE);cos(azStartE);0];

%The initial target state with GLOBAL position and LOCAL velocity is thus
xInit=[rGlobalCart;vLocalInit;0.004];

times = linspace(0,1200,100);

%The continuous-time drift function for a 3D, flat-Earth dynamic model is
aDrift=@(x,t)aCircTurn(x);

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
[xList] = StateProp3D(xInit, times, aDrift);

% % Plots
figure
hold on, grid on, box on
plot3(xList(1,:), xList(2,:), xList(3,:))

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
visible = coast * mean(xList(1:3,:),2) >= 0;  % Set visibility condition (z >= 0 is near side)

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
