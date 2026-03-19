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

% Doug Carlson (doug.o.carlson@gmail.com), 2026-03-18 20:54

%%

fc = 10*GHz;
BW = 10*MHz;

Pt = 100*kW;
PRF = 1*kHz;

tauP = 20/BW;

% X-band LFM
tx = Radar.RadarTransmitter(fc, Pt, tauP, PRF, ...
    'bandwidth', BW, ...
    'waveform', 'lfm');

tx   % calls disp()

Fs = 2*BW;
Ts = 1/Fs;
nFFT = 256;

[s, t]  = tx.generatePulse(Fs);
[Pout_dB, BW] = tx.spectrum(Fs, nFFT);

h       = tx.matchedFilter(Fs);
rc      = conv(s, h, 'full');    % range-compressed pulse

% % Barker-13
% b13 = pi * [0 0 0 0 0 1 1 0 0 1 0 1 0]';
% txB = Radar.RadarTransmitter(3e9, 10e3, 26e-6, 500, ...
%     'waveform', 'barker', 'phaseCode', b13);
[chi, tau_ax, fd_ax] = tx.ambiguityFunction(Fs, nFFT);
[FD, TAU] = meshgrid(fd_ax, tau_ax);

%
figure
surf(TAU,FD,chi)
shading flat
colorbar