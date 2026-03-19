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
% Doug Carlson (doug.o.carlson@gmail.com), 2024-10-27 21:09

%%

% Parameters for the clock display
radius = 1;            % Radius of the clock face
center = [0, 0];       % Center of the clock
numSeconds = 60;       % Number of seconds in a full rotation

% Create figure for the clock
figure;
hold on;
axis equal;
axis([-1.5, 1.5, -1.5, 1.5]);
title('Clock Second Hand Simulation');

% Draw clock circle
theta = linspace(0, 2*pi, 100);
plot(radius*cos(theta), radius*sin(theta), 'k'); % Clock face border

% Animate the second hand
for second = 0:numSeconds-1
    % Calculate the angle of the second hand (clockwise direction)
    angle = -2 * pi * second / numSeconds;
    
    % Calculate the end point of the second hand
    x = radius * cos(angle);
    y = radius * sin(angle);
    
    % Plot the second hand
    secondHand = plot([center(1), x], [center(2), y], 'r', 'LineWidth', 2);
    
    % Pause for 1 second
    pause(1);
    
    % Delete the old second hand
    delete(secondHand);
end
