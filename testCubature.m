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

% Doug Carlson (doug.o.carlson@gmail.com), 2026-04-14 20:37

%%

zDim = 1;

%Cubature points for measurement conversion.
[xi,w] = fifthOrderCubPoints(zDim);


Ns = 2^20;
m = 0;
S = 1;
z = m + S*randn(zDim,Ns);
mean(z)
std(z)

% Nonlinear transform
h = @(x) x.^2/2;

mean(h(z))
var(h(z))

xit = h(transformCubPoints(xi,m,S));
[mu, P] = calcMixtureMoments(xit,w)

figure
tiledlayout(2,1)
nexttile
histogram(z,Normalization="pdf")
nexttile
hold on
histogram(h(z), Normalization="pdf")

x = 1e-3:1e-2:20;
y = exp(-x)./sqrt(pi*x);
plot(x,y)
xlim([0,2])

