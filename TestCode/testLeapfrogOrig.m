clc
clear
close all

cd(fileparts(which(mfilename)))

addpath C:/Users/dougo/OneDrive/Documents/MATLAB/
Constants

dbstop if error

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-10-25 17:36

%% Rotator Target Intercept Using Leapfrog Method
clc
% Number of samples
Ns = 10;

% Propagator
F = @(t) kron([1,t;0,1],eye(3));
f = @(t,s) F(t)*s;

% Azimuth
az = @(s) mod(atan2d(s(1,:),s(2,:)), 360);

% Rotator
Ts = 10;
pos_rot_enu = [0;0;0];

% Target
pos_tgt_enu = [-10;0;0]*km;
vel_tgt_enu = 125*[1;1;0]/sqrt(2);
s_enu = [pos_tgt_enu; vel_tgt_enu] * ones(1, Ns);
t = zeros(1, Ns);

az0 = 0;
for k = 1:Ns
    iter = 0;
    while iter < 5
        az1 = az(s_enu(:,k));
        dt = Ts*(az1-az0)/360;
        s_enu(:,k) = f(dt,s_enu(:,k));
        t(k) = t(k) + dt;
        az0 = az1;
        iter = iter+1;
    end
    if k < Ns
        s_enu(:,k+1) = s_enu(:,k);
        t(k+1) = t(k);
        az0 = az0-360;
    end
end
t
s_enu

%% Plots
% Range ring
r0 = 15;
th = (0:360)*deg;

Ms = 4000*pi;

figure(Position=[480 130 560 420])
hold on, grid on, box on
plot(s_enu(1,:)/km, s_enu(2,:)/km, 'ko')
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

    s = f(Ns*Ts/Ms,s);
end

