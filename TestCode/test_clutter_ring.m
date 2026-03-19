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
% Doug Carlson (doug.o.carlson@gmail.com), 2025-06-09 20:15

%%
pos_plat_llh = [30*deg;-82*deg;3*km];
% pos_plat_llh = [ 0*deg; 0*deg;3*km];
pos_plat_ecf = ellips2Cart(pos_plat_llh);
CE = getENUAxes(pos_plat_llh, false);

% Slant range
rs = 150*km;

% Clutter azimuth
Nc = 36;
phc = 2*pi*(0:Nc-1)/Nc;
% phc = 0
OneC = ones(size(phc));

% Guess
thN = acos(pos_plat_llh(3)/rs) * OneC;

a = Constants.WGS84SemiMajorAxis;
b = Constants.WGS84SemiMinorAxis;
A = diag([a,a,b].^2)\eye(3);

max_iter = 10;
iter = 0;
abs_del = inf;
tol = 1e-6;
while (iter<max_iter) && (abs_del > tol)
    P = pos_plat_ecf + rs*CE*[sin(phc).*sin(thN);cos(phc).*sin(thN);-OneC.*cos(thN)];
    dP = rs*CE*[sin(phc).*cos(thN);cos(phc).*cos(thN);OneC.*sin(thN)];

    fN = sum(P .* (A*P), 1) - 1;
    dfN = 2*sum(dP .* (A*P), 1);

    del = fN./dfN;
    abs_del = max(abs(del));
    thN = thN - del;
    iter = iter+1;
end

iter
abs_del
pos_patch_ecf = pos_plat_ecf + rs*CE*[sin(phc).*sin(thN);cos(phc).*sin(thN);-OneC.*cos(thN)];
pos_patch_llh = Cart2Ellipse(pos_patch_ecf);

% fig = figure;
clf
gx = geoaxes();
hold on

geoplot(gx,pos_plat_llh(1)/deg,pos_plat_llh(2)/deg,'o')
geoplot(gx,pos_patch_llh(1,:)/deg,pos_patch_llh(2,:)/deg,'k.')
geobasemap(gx,'colorterrain')



