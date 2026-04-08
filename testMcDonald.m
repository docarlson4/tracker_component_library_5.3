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

% Doug Carlson (doug.o.carlson@gmail.com), 2026-04-04 15:50

%% References
% [1] McDonald, M. K., & Cerutti-Maori, D. (2016). Coherent radar
% processing in sea clutter environments, part 1: modelling and partially
% adaptive STAP performance. IEEE Transactions on Aerospace and Electronic
% Systems, 52(4), 1797-1817.

%% Constants
% Create a sphere
a = Constants.WGS84SemiMajorAxis;
b = Constants.WGS84SemiMinorAxis;

%% Platform Info

% Spectrum
fc = 9.09*GHz; % 9.09 and 9.50 GHz
lamc = c0/fc;
BW = 20*MHz;
PRF = 3*kHz;
CPI = 42.6*ms;
Np = ceil(CPI*PRF);

% Geometry
hp = 2.5*km;
vp = 100;
psiG = 20*deg;

% Slant range
rs = -a*sin(psiG) + sqrt((a*sin(psiG))^2 + 2*a*hp + hp^2);
% Ground range
phi = acos((a^2 + (a+hp)^2 - rs^2)/(2*(a+hp)*a));
rho = a*phi;

% Pass 1 - downwind (radar look direcion in downwind dir)
% Pass 2 - upwind (radar look direcion into wind)

% Antenna
Nc = 3;
antLen = 80*cm;
antPol = 'VV';

% Scene center
pos_scn_geo = [53.98*deg; 7.91*deg; 0];
pos_scn_ecf = ellips2Cart(pos_scn_geo);
CE = getENUAxes(pos_scn_geo);

% Sea conditions 
% - Sea conditions were moderate, with a reported 
%   - swell height of 0.9–1.5 m
%   - wind sea wave height of 0.4–0.5 m.


