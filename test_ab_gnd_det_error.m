clc
clear
close all

cd(fileparts(which(mfilename)))
addpath c:\Users\dougo\OneDrive\Documents\MATLAB\
Constants

dbstop if error

%% Airborne (AB) Monostatic Ground Target Detection Error
num_samples = 1000;

%% Aircraft Position
pos_orig_llh = [34*deg; -84*deg];

alt_plat = 6*km;
pos_plat_llh = [pos_orig_llh; alt_plat];

%% Radar Parameters
f = 10*GHz;
lam = c0/f;
len_arr = 1;

SNR = 10;
snr = invdB(SNR);

B = 40*MHz;
rng_res = c0/(2*B);
sig_rng = rng_res/sqrt(2*snr);

az_bw = sqrt(pi/4) * lam/len_arr;
az_res = az_bw;
sig_az = az_res/sqrt(2*snr);

%% Elevation (wrt Nadir) Error

% DEM error
sig_hgt = 1;

% Target slant range: rs
rs = 20*km + sig_rng*randn(1,num_samples);

% Target bearing (true-North)
az = 00*deg + sig_az*randn(1,num_samples);

% Target elevation
a2 = alt_plat + RE + sig_hgt*randn(1,num_samples);
a1 = RE + sig_hgt*randn(1,num_samples);
num = a2.^2 - a1.^2 + rs.^2;
den = 2*a2.*rs;
cos_el = num./den;
el = acos(cos_el);

% Target ground range
del = asin(sin(el).*rs./a1);
rho = a1.*del;

% Final latitude/longitude position
[pos_tgt_llh] = destination(pos_orig_llh, az, del);

sig_lat_lon = std(pos_tgt_llh,[],2)
sig_E_N = [
    sig_lat_lon(2)*sin(mean(pos_tgt_llh(1,:)));
    sig_lat_lon(1)] * RE

points = [pos_tgt_llh;zeros(1,num_samples)];
cartPoints = ellips2Cart(points,RE,0);
[xENU,M]=ECEF2ENU(mean(points,2),cartPoints,RE,0);
std(xENU,[],2)

%% Helper Functions
function [ll] = destination(ll0, az, del)

ll(1,:) = asin(sin(ll0(1))*cos(del) + cos(ll0(1))*sin(del).*cos(az));
ll(2,:) = ll0(2) + atan2( ...
    sin(az).*sin(del)*cos(ll0(1)), cos(del) - sin(ll0(1))*sin(ll(1)) );

end





