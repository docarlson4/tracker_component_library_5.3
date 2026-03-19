clc
clear
close all

cd(fileparts(which(mfilename)))
startup_tcl

dbstop if error

MyConstants

% Developed in Matlab 9.12.0.2327980 (R2022a) Update 7 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2023-12-03 17:04

%% 

% Number of trials
num_trials = 1e5;

% Number of resolution cells: n
n = 1000;

% Probability of false alarm
Pfa = 0.93

% Expected number of false alarms
lamV = n*(1 - (1-Pfa)^(1/n))

% Bernoulli trials
p = lamV/n % that cell is occcupied
r = rand(n, num_trials);
b = sum(r <= p);

figure
subplot(211)
histogram(b)
a=axis;

% Pfa = 1 - (1-p)^n

Pfa_experiment = sum(b>0)/num_trials

pois = PoissonD;
n_pois = pois.rand([num_trials,1], lamV);
% n_pois = poissrnd(lamV,num_trials);
subplot(212)
histogram(n_pois)
axis(a)

