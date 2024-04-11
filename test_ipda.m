clc
clearvars
close all

cd(fileparts(which(mfilename)))

startup_tcl

dbstop if error

MyConstants

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-04-10 19:20

%% Monte Carlo Samples
Nmc = 1;

%% Measurement Model
mmt_space = [
    0, 0
    1e3, 400];
mmt_vol = prod(diff(mmt_space)); % m^2

lambda = 1e-5; % 1/scan/m^2
lambdaV = lambda * mmt_vol;

PD = 0.9;
sig_xy = 5;
SR = diag([sig_xy,sig_xy]);
R = SR*SR;

H = kron([1,0],eye(2));

%% Motion Model
T = 1; % s
numSamples = 21;
t = (0:numSamples-1)*T;

F1 = [
    1,T
    0,1];
F = kron(F1,eye(2));

Q1 = [
    T^4/4,T^3/2
    T^3/2,T^2];
q = 0.75;
Q = kron(Q1,eye(2));
SQ = cholSemiDef(Q);

% Markov Chain 1
p11 = 0.98; p21 = 0;
p12 = 0.02; p22 = 1;

% Track existence probability
PIsRealInit = 0.2;

% Gate Probability
PG = 0.98;
gammaVal = ChiSquareD.invCDF(PG,2);


%% Truth Data
xDim = 4;
s0 = [130;200;35;0] * ones(1,numSamples);
for curScan = 2:numSamples
    s0(:,curScan) = F*s0(:,curScan-1);
end

%% Measurements
zMeasCart = cell(numSamples,Nmc);
for kMC = 1:Nmc

    for curScan = 1:numSamples

        numFalse = PoissonD.rand(1,lambdaV);
        numTargets = 1;
        if(curScan <= numSamples)
            %Determine which, if any, targets should be detected.
            isDet = rand(numTargets,1) < PD;
        else
            %We are after the end; the tracks are now gone.
            isDet = [0;0];
        end

        %Allocate space for the detections.
        numMeas = numFalse+sum(isDet);
        zCartCur = zeros(2,numMeas);
        curDet = 1;

        if(numMeas>0)
            %Generate the detection from the targets, if any.
            for curTar = 1:numTargets
                if(isDet(curTar))

                    zCartCur(:,curDet) = H*s0(:,curScan) + SR*randn(2,1);

                    curDet = curDet+1;
                end
            end

            %Generate the false alarm detections, if any.
            xClutBounds = mmt_space(:,1);%[minRange;maxRange];
            yClutBounds = mmt_space(:,2);%[-pi;pi];
            for curFalse = 1:numFalse
                x = UniformD.rand(1,xClutBounds);
                y = UniformD.rand(1,yClutBounds);

                zCartCur(:,curDet) = [x;y];
                curDet = curDet+1;
            end

            %We will now convert the measurements into Cartesian coordinates as
            %we are using a converted-measurement filter.
            zMeasCart{curScan, kMC} = zCartCur;
        end

    end

end

%% Filter
for kMC = 1:Nmc
%We will save the value of each track at each time. Thus, a structure
%holds the track states at each time. The first field in state_st are the
%actual values, the fourth field in state_st is an ID so that we can
%associate the tracks across time.

state_st = struct('x',[],'S',[],'r',[],'ID',[]);
state_st = repmat(state_st, [numSamples,1]);

track_st = struct('x',[],'S',[],'r',[],'ID',[],'scan_num',[],'num_hits',[]);
   
clear two_point_init
for curScan = 1:numSamples
    tCur = (curScan-1)*T;
    zCur = zMeasCart{curScan};
    SRCur = SR;
    numMeas = size(zCur,2);
    % case "Two-Point"
    [xNew, SNew] = two_point_init(zCur, SRCur, tCur, [25,50], xDim);

    numStates = size(xNew, 2);

    %Initialization existence probabilities
    rNew = PIsRealInit*ones(1,numStates);

    %Next, if there are any existing states, we want to predict them to the
    %current time step and update them with the measurements
    if ( curScan > 1 ) && ~isempty( state_st(curScan-1).x )
        x = state_st(curScan-1).x;
        S = state_st(curScan-1).S;
        r = state_st(curScan-1).r;
        xID = state_st(curScan-1).ID;
        numTargetsCur = size(x,2);
        
        %Predict the tracks to the current time.
        for curTar = 1:numTargetsCur
            [x(:,curTar),S(:,:,curTar)] = sqrtDiscKalPred( ...
                x(:,curTar),S(:,:,curTar),F,SQ);
        end
        %Update the target existence probabilities with the Markov 
        %switching model.
        r = p11*r;

        %The inclusion of r takes into account the track existence
        %probabilities.
        measJacobDets = [];
        [A,xHyp,PHyp] = makeStandardCartOnlyLRMatHyps(x,S,zCur,SRCur,[], ...
            PD,lambda,r,gammaVal,measJacobDets);

        disp('')
    else
        state_st(curScan).x = xNew;
%         state_st(curScan).ID = IDNew;
        state_st(curScan).S = SNew;
        state_st(curScan).r = rNew;        
    end    

    %% Populate track_st
    track_st = get_tracks(curScan, state_st(curScan), track_st);

end

end


%% Plots

figure('Position', [300 340 1000 420])

hold on, grid on, box on

for kS = 1:numSamples
    if ~isempty(zMeasCart{kS,1})
        scatter(zMeasCart{kS,1}(1,:),zMeasCart{kS,1}(2,:),'.')
    end
end
plot(s0(1,:),s0(2,:),'ko')

axis equal
axis(mmt_space(:)')











