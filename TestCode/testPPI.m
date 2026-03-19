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
% Doug Carlson (doug.o.carlson@gmail.com), 2024-10-27 21:07

%%

% Parameters for PPI display
numRadialBins = 200;  % Number of range bins
numAzimuthBins = 360; % Number of azimuth bins (0-360 degrees)
maxRange = 100;       % Maximum range in arbitrary units

% Generate synthetic radar data (example, could be replaced by actual data)
radarData = rand(numRadialBins, numAzimuthBins); % Simulated radar returns

% Set up the polar coordinate grid
theta = linspace(0, 2*pi, numAzimuthBins);  % Azimuth angles in radians
r = linspace(0, maxRange, numRadialBins);   % Range values

% Create a meshgrid for the polar coordinates
[Theta, R] = meshgrid(theta, r);

% Convert polar coordinates to Cartesian for plotting
[X, Y] = pol2cart(Theta, R);

% Plot the PPI display
figure;
pcolor(X, Y, radarData);
shading interp;               % Smooth shading for a continuous look
colormap(jet);                % Set colormap (adjust as desired)
colorbar;                     % Add colorbar to indicate intensity
axis equal;                   % Equal scaling on both axes
title('PPI Display');
xlabel('X Range');
ylabel('Y Range');
