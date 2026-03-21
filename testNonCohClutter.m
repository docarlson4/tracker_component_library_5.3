clc
clearvars
close all

cd(fileparts(which(mfilename)))

restoredefaultpath

addpath C:/Users/dougo/OneDrive/Documents/MATLAB/
Constants

addpath(genpath("C:/Users/dougo/source/repos/tracking/tracker_component_library_5.3/"))
Constants;

dbstop if error

% Doug Carlson (doug.o.carlson@gmail.com), 2026-03-20 13:47

%% Spectrum
fc = 10*GHz;
lamc = c0/fc;
B = 100*MHz;
rho = c0/B;

%% Radar Range
R = 50*km;

%% Radar Antenna
Lant = 1;
the3dB = sqrt(pi/4)*lamc/Lant;

%% Clutter Patch
A = 0.75 * rho * R * the3dB;
nu_obj = SeaClutter.ShapeParameter("ResCellArea", A);

sc_obj = SeaClutter.SeaClutter('Shape', nu_obj.Shape);
sc_obj.GenClutter
sc_obj.PlotTexture("Image")

nc_obj = SeaClutter.NonCoherentClutter(CorrCoeffRng=@(x,tau) exp(-abs(x)/tau))
