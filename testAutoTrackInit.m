clc
clearvars
close all

startup_tcl

MyConstants

dbstop if error

%% Radar Parameters
radar_obj = RadarReceiver( ...
    "Frametime", 2, ...
    "Bandwidth", 10*MHz, ...
    "NumPulse", 20, ...
    "SNRdB", 10, ...
    "PRF", 1/(300*us));

% radar_obj.PlotUnfoldingProb

%% Model Setup
disp('Setting parameters and creating the simulated trajectories.')

zDim = 2;

rng(0)
% rng("shuffle")

PD = 0.99998;%Detection probability --same for all targets.

% Probability of false target
PFT = 0.99995;

%The viewing region range. This is important for dealing with the clutter
%model.
RangeMin = 1e3;
RangeMax = radar_obj.RangeAmb;
%The angular region is two pi (all around).
mmt_space = [
    RangeMin 0 -radar_obj.RangeRateAmb/2
    RangeMax  2*pi  radar_obj.RangeRateAmb/2];
% Measurement volume
mmt_vol = prod(diff(mmt_space));

% Elliptical gate probability
PG = 0.95;
gammaVal = ChiSquareD.invCDF(PG,zDim);

%% Generate Target Truth Data

Ts = radar_obj.FrameTime;
NumSamples = 60;
SpeedMin = 205*mph;
SpeedMax = 600*mph;
NumTgts = 5;

% [ZCartTrue, ZPolTrue, Ts, num_tgt_truth, max_acc] = ...
%     genSwarm(Ts, NumSamples, RangeMin, RangeMax, SpeedMin, SpeedMax, NumTgts);
% numSamples = length(ZCartTrue);

acs_obj = AircraftSwarm( ...
    "Type", "FixedWing", ...
    "NumAircraft", NumTgts, ...
    "NumSamples", NumSamples, ...
    "MeanLifetime", 75, ...
    "MeanFirstDet", 10, ...
    "FrameTime", radar_obj.FrameTime, ...
    "RangeMax", RangeMax, ...
    "RangeMin", RangeMin, ...
    "SpeedMax", SpeedMax, ...
    "SpeedMin", SpeedMin);
ZCartTrue = acs_obj.ZCart;
ZPolTrue = acs_obj.ZPol;

%% Generate Measurements
disp('Generating measurements.')

%lambda times the "volume" in the receiver's polar coordinate system needed
%for the Poisson clutter model
lambdaV = -log( 1 - PFT );

%Assumed standard deviations of the measurement noise components.
sigmaR = radar_obj.RangeUnc;
sigmaAz = radar_obj.AzimuthUnc;
sigmaRr = radar_obj.RangeRateUnc;
%Square root measurement covariance matrix; assume no correlation.
SR = diag([sigmaR,sigmaAz,sigmaRr]);

[zMeasCart, SRMeasCart, zMeasPol] = genMmts( ...
    ZCartTrue, ZPolTrue, PD, lambdaV, mmt_space, radar_obj);
%     ZCartTrue, ZPolTrue, PD, lambdaV, SR, mmt_space, Ts, radar_obj.RangeRateAmb);

%% Use of Range-Rate Measurements in Automatic Track Formation

figure(Position=[207.4000  152.2000  840.8000  609.6000])
hold on, grid on, box on

dkHist = [];
for kS = 2:NumSamples

    % Filter stationary clutter (SCF)
    r1 = zMeasPol{kS-1}(1,:);
    r2 = zMeasPol{kS}(1,:)';
    a1 = zMeasPol{kS-1}(2,:);
    a2 = zMeasPol{kS}(2,:)';
    [i2,i1] = find(sqrt(r2.*r2-2*r2.*r1.*cos(a2-a1)+r1.*r1)<50);
    i1u = unique(i1)';
    i2u = unique(i2)';
    r1(i1u) = nan;
    lgcl1 = ~isnan(r1);
    r2(i2u) = nan;
    lgcl2 = ~isnan(r2);
    
    % Gather remaining from SCF
    r1 = zMeasPol{kS-1}(1,lgcl1);
    rr1 = zMeasPol{kS-1}(3,lgcl1);
    t1 = zMeasPol{kS-1}(4,lgcl1);
    n1 = length(t1);

    r2 = zMeasPol{kS}(1,lgcl2)';
    rr2 = zMeasPol{kS}(3,lgcl2)';
    t2 = zMeasPol{kS}(4,lgcl2)';
    n2 = length(t2);

    % Estimate wrapping factors
    a1 = round(( (r2-r1)./(t2-t1) - rr1 )/radar_obj.RangeRateAmb );
    a2 = round(( (r2-r1)./(t2-t1) - rr2 )/radar_obj.RangeRateAmb );

    % Perform least-squares estimate of range and range-rate
    Tk = t2-t1; Tk = Tk(:)';
    rm1 = ones(n2,1) * r1; rm1 = rm1(:)';
    rrm1 = ones(n2,1) * rr1 + a1 * radar_obj.RangeRateAmb; rrm1 = rrm1(:)';
    rm2 = r2 * ones(1,n1); rm2 = rm2(:)';
    rrm2 = rr2 * ones(1,n1) + a2 * radar_obj.RangeRateAmb; rrm2 = rrm2(:)';
    Yk = [rm1;rrm1;rm2;rrm2];

    SRk = kron(eye(2),SR([1,3],[1,3]));
    Rk = SRk*SRk';

    Xk = zeros(2,n1*n2);
    dk = zeros(1,n1*n2);
    for k = 1:n1*n2
        Hk = [[1 -Tk(k); 0 1];eye(2)];
        Jk = Hk' * (Rk\Hk);
        Xk(:,k) = (Jk\Hk') * (Rk\Yk(:,k));
        zk = SRk\(Yk(:,k) - Hk*Xk(:,k));
        dk(k) = zk'*zk;
    end

    % Filter according to Chi-sqr statistic
    lgclKeep = dk < 1*gammaVal;
    idxKeep = find(lgclKeep);
    [i2,i1] = ind2sub([n2,n1],idxKeep);
    i2u = unique(i2);
    dkHist = [dkHist, dk(lgclKeep)]; %#ok<AGROW> 

    % Scatter plots to visually gauge performance
    scatter(zMeasCart{kS}(1,:)/km, zMeasCart{kS}(2,:)/km, '.')
    scatter(zMeasCart{kS}(1,i2u)/km, zMeasCart{kS}(2,i2u)/km, 'o')
    axis equal
    axis([-1,1,-1,1]*50)
    drawnow
end

%%
figure
hold on, grid on, box on
histogram(dkHist, 'Normalization','pdf')
histogram(ChiSquareD.rand([length(dkHist),1],2), 'Normalization','pdf')
disp('')

return
%%
% MATLAB Code to Wrap True Range-Rate Values according to Ambiguous Range-Rate

% Define the ambiguous range-rate (ΔRR)
drr = 100; % Example value, you can change this to your specific drr

% Define the set of true range-rate values (rrt)
rrt = [-150, -90, -60, 0, 45, 120, 200]; % Example values, you can change these

% Compute the wrapped range-rate values
wrapped_rrt = mod(rrt + drr/2, drr) - drr/2;

% Display the results
fprintf('True Range-Rate Values (rrt):\n');
disp(rrt);
fprintf('Wrapped Range-Rate Values (wrapped_rrt):\n');
disp(wrapped_rrt);





