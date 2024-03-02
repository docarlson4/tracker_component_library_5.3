clc
clear
close all

cd(fileparts(which(mfilename)))
addpath Utilities\

dbstop if error

Constants

% Developed in Matlab 9.12.0.2327980 (R2022a) Update 7 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-03-02 13:55

%%
Ts = 12;
Ns = 100;
Rmin = 1*km;
Rmax = 40*km;
Vmin = 25*mph;
Vmax = 600*mph;
lam = 15;
[ZCartTrue, ZMeasTrue, Ts, num_tgt] = genSwarm(Ts, Ns, Rmin, Rmax, Vmin, Vmax, lam);

figure
hold on, grid on, box on
for kS = 1:Ns
    if ~isempty(ZCartTrue{kS})
        scatter(ZCartTrue{kS}(2,:)/km, ZCartTrue{kS}(1,:)/km)
    end
    axis equal
    axis([-1,1,-1,1]*2*Rmax/km)
    drawnow
end
