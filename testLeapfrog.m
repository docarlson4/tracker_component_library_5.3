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
% Doug Carlson (doug.o.carlson@gmail.com), 2024-10-25 17:36

%% Rotator Target Intercept Using Leapfrog Method
clc
% Number of samples
Ns = 10;

% Propagator
F = @(t) kron([1,t;0,1],eye(3));
fProp = @(t,s) F(t)*s;

% Azimuth
az = @(s) mod(atan2d(s(1,:),s(2,:)), 360);

% Rotator
Ts = 10;
pos_rot_enu = [0;0;0];

% Target
pos_tgt_enu = [-10;0;0]*km;
vel_tgt_enu = 125*[1;1;0]/sqrt(2);
sInitENU = [pos_tgt_enu; vel_tgt_enu];

% Target-Beam intercept
[sFinalENU, t] = LeapFrogENU(sInitENU, Ns, Ts, fProp);

%% Plots
% Range ring
r0 = 15;
th = (0:360)*deg;

Ms = 4000*pi;

figure(Position=[480 130 560 420])
hold on, grid on, box on
plot(sFinalENU(1,:)/km, sFinalENU(2,:)/km, 'ko')
plot(r0*sin(-th), r0*cos(th))
axis equal
axis([-10,10,-10,10]*1.5)
ax = gca;

s0 = [pos_tgt_enu; vel_tgt_enu];
s = [pos_tgt_enu; vel_tgt_enu];
for k = 0:Ms-1

    plt0 = plot(s(1)/km, s(2)/km, 'r.');
    plt1 = plot(ax, ...
        [0,r0*sin(2*pi*k*Ns/Ms)], ...
        [0,r0*cos(2*pi*k*Ns/Ms)],'k');
%     pause(0.1)
    drawnow
    delete([plt0,plt1]);

    s = fProp(Ns*Ts/Ms,s);
end

