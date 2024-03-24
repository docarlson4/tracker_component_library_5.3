clc
clear
close all

cd(fileparts(which(mfilename)))

dbstop if error

MyConstants

% Developed in Matlab 9.12.0.2327980 (R2022a) Update 7 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-03-13 14:13

%% Model Setup
disp('Setting parameters and creating the simulated trajectories.') 

% rng("shuffle")
rng(0)

PD = 0.98;%Detection probability --same for all targets.

% Probability of false target
PFT = 0.99995;

%lambda times the "volume" in the receiver's polar coordinate system needed
%for the Poisson clutter model
lambdaV = -log( 1 - PFT );

%The viewing region range. This is important for dealing with the clutter
%model.
minRange = 1e3;
maxRange = 40e3;
%The angular region is two pi (all around).

%Cubature points for measurement conversion.
zDim = 2;
[xi,w] = fifthOrderCubPoints(zDim);


%% Generate Target Truth Data

Ts = 12;
Ns = 100;
Rmin = minRange;
Rmax = maxRange;
Vmin = 25*mph;
Vmax = 600*mph;
lam = 3;

[ZCartTrue, ZPolTrue, Ts, num_tgt_truth] = ...
    genSwarm(Ts, Ns, Rmin, Rmax, Vmin, Vmax, lam);
numSamples = length(ZCartTrue);
numSamplesTrack = numSamples;

%% Generate Measurements
disp('Generating measurements.') 

%Assumed standard deviations of the measurement noise components.
sigmaR = 10;
sigmaAz = 0.1*deg;
%Square root measurement covariance matrix; assume no correlation.
SR = diag([sigmaR,sigmaAz]);

%Generate measurements and false alarms for each scan.
zMeasCart = cell(numSamples,1);
SRMeasCart = cell(numSamples,1);
zMeasJacobDet = cell(numSamples,1);
for curScan = 1:numSamples
    %Determine the number of false alarms to generate.
    numFalse = PoissonD.rand(1,lambdaV);
    numTargets = size(ZCartTrue{curScan},2);
    if(curScan <= numSamplesTrack)
        %Determine which, if any, targets should be detected.
        isDet = rand(numTargets,1) < PD;
    else
        %We are after the end; the tracks are now gone.
        isDet = [0;0];
    end

    %Allocate space for the detections.
    numMeas = numFalse+sum(isDet);
    zCur = zeros(2,numMeas);
    curDet = 1;

    if(numMeas>0)
        %Generate the detection from the targets, if any.
        for curTar = 1:numTargets
            if(isDet(curTar)) && ( ZPolTrue{curScan}(1,curTar) < maxRange)
                zCur(:,curDet) = ZPolTrue{curScan}(:,curTar) ...
                    + SR*randn(2,1);
                curDet = curDet+1;
            end
        end

        %Generate the false alarm detections, if any. 
        rClutBounds = [minRange;maxRange];
        azBounds = [-pi;pi];
        for curFalse = 1:numFalse
            r = UniformD.rand(1,rClutBounds);
            az = UniformD.rand(1,azBounds);

            zCur(:,curDet) = [r;az];
            curDet = curDet+1;
        end
            
        %We will now convert the measurements into Cartesian coordinates as
        %we are using a converted-measurement filter.
        [zMeasCart{curScan},RMeasCart] = pol2CartCubature(zCur,SR,0,true, ...
            [],[],[],xi,w);
        
        %Take the lower-triangular square root of the covariance matrices.
        measDetCur = zeros(numMeas,1);
        for curMeas = 1:numMeas
            RMeasCart(:,:,curMeas) = chol(RMeasCart(:,:,curMeas),'lower');
            measDetCur(curMeas) ...
                = det(calcPolarConvJacob(zCur(:,curMeas),0,true));
        end
        SRMeasCart{curScan} = RMeasCart;
        zMeasJacobDet{curScan} = measDetCur;
    end
end


