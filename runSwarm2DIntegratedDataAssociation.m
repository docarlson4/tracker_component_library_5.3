%%DEMO2DINTEGRATEDDATAASSOCIATION Demonstrate how the Global Nearest
%                       Neighbor (GNN) Joint Integrated Probabilistic Data
%                       Association Filter (JIPDAF) can be used to track
%                       targets performing automatic track initiation and
%                       termination. This is a simple scenario where no
%                       gating/ clustering is performed, so approximate
%                       data association probabilities must be used for the
%                       problem to be computationally practicable.
%
%This is a simple two-dimensional (x,y) simulation scenario involving two
%maneuvering aircraft that cross, come within close range of each other and
%then separate. The scenario involves the trajectory used in the air
%traffic control (ATC) scenario discussed in Chapter 11.7.4 of [1] and a
%shifted, mirrored version of it. As measurement filtering is not the focus
%of this problem, a converted-measurement Kalman filter with a first-order
%Gauss-Markov (integrated Ornstein-Uhlenbeck) dynamic model is used with
%polar-Cartesian converted measurements.
%
%Such a simulation scenario might be representative of what one gets from a
%rotating radar when all contacts from a single rotation are collected into
%a scan rather than performing updates for every single dwell. In such a
%scenario, not all of the measurements will be taken at the same time,
%meaning that different targets will be predicted to different times based
%on where they are. To simplify the presentation here, all measurements are
%taken to be at the same time.
%
%As approximate target-measurement association probabilities are used here,
%it was observed that the algorithm is sufficiently fast without performing
%gating and then jointly processing groups of targets that gate with common
%measurements together. That is, all targets are processed jointly in this
%implementation, which also simplifies it. If gating and grouping targets
%are neded to reduce computational complexity in other applications (such
%as if one wishes to use exact target-measurement association
%probabilities), then one can use the DisjointSetM class to cluster the
%targets and measurements into groups. Such clustering, without discussion
%of data structures to perform it, is mentioned in the original paper on
%the JIPDA in [2].
%
%The scenario is run with false alarms and missed detections. To illustrate
%timely termination of tracks, an extra 50 time-steps of simulation are
%performed after the trajectories end to demonstrate that they are
%terminated in a timely manner.
%
%This function plots all of the tracks found by the tracking algorithm in
%green. Tracks are only initially displayed once their existence
%probability exceeds 95% (the PDisp parameter). After that, they continue
%to be displayed until they are terminated. Tracks are terminated if their
%probability of existence is less than 0.01%.
%
%Note that performance could be improved if range-rate (Doppler)
%information were simulated and used. Specifically, even if not used in the
%measurement update step, range-rate information could improve the
%likelihoods used for target-measurement association, reducing the
%occurrence of false tracks.
%
%REFERENCES:
%[1] Y. Bar-Shalom, X. R. Li, and T. Kiruabarajan, Estimation with
%    Applications to Tracking and Navigation. New York: Wiley Interscience,
%    2001.
%[2] D. Musicki and R. Evans, "Joint integrated probabilistic data
%    association: JIPDA," IEEE Transactions on Aerospace and Electronic
%    Systems, vol. 40, no. 3, pp. 1093-1099, Jul. 2004.
%
%July 2016 David F. Crouse, Naval Research Laboratory, Washington D.C.
%(UNCLASSIFIED) DISTRIBUTION STATEMENT A. Approved for public release.

clc
clearvars
% close all

startup_tcl

MyConstants

dbstop if error

km = 1e3;

disp('An example of a simple Kalman filter with Gauss-Markov model, data association, and track initiation/ termination in 2D') 

%% Model Setup
disp('Setting parameters and creating the simulated trajectories.') 

zDim = 2;

rng(0)
% rng("shuffle")

PD = 0.98;%Detection probability --same for all targets.

% Probability of false target
PFT = 0.95;

%lambda times the "volume" in the receiver's polar coordinate system needed
%for the Poisson clutter model
lambdaV = -log( 1 - PFT );
% lambdaV = 2;

%The viewing region range. This is important for dealing with the clutter
%model.
minRange = 1e3;
maxRange = 40e3;
%The angular region is two pi (all around).
mmt_space = [minRange -pi; maxRange, pi];
% Measurement volume
mmt_vol = prod(diff(mmt_space));

%Assumed clutter parameter in terms of false alarms per meter-radian
%(polar coordinates). False alarms are generated in the measurement
%domain of the radar. The tracker can still work if the value for lambda
%does not perfectly match the value simulated/ the real value. However, one
%should note that lambda in the tracker cannot be arbitrarily large lest
%the tracker eventually ignore all measurements (false alarms become always
%more likely than true tracks).
lambda0 = lambdaV/mmt_vol;

% Clutter estimation class (lambda equivalent)
Lx = 2*maxRange; lx = 1*km;
Ly = 2*maxRange; ly = 1*km;
cm_obj = Tracker.ClutterEstimation("Type","Classical", ...
    "AveragingLength", 25, "MmtRegion", [-Lx/2, -Ly/2; Lx/2, Ly/2], ...
    "CellSize", [lx,ly]);

%The AlgSel parameters are inputs to the singleScanUpdate function.
%Use the GNN-JIPDAF with approximate association probabilities.
algSel1 = 3;
algSel2 = 0;
param3 = [];

%The transition probability when predicting forward the target states that
%the track did not end during the transition. This is a design parameter. 
PStayAlive = 0.98;

%This is the initial probability assigned that the detection might be a
%real, new track. This value is a design parameter.
PIsRealInit = 0.1;

%The probability of existence below which a track is terminated.
PTerminate = 1e-4;

% Elliptical gate probability
PG = 0.95;
gammaVal = ChiSquareD.invCDF(PG,zDim);

%% Generate Target Truth Data

Ts = 12;
Ns = 100;
Rmin = minRange;
Rmax = maxRange;
Vmin = 25*mph;
Vmax = 600*mph;
vLims = [Vmin, Vmax];
avgNumTgts = 3;

[ZCartTrue, ZPolTrue, Ts, num_tgt_truth, max_acc] = ...
    genSwarm(Ts, Ns, Rmin, Rmax, Vmin, Vmax, avgNumTgts);
numSamples = length(ZCartTrue);
% maxAccel = max(max_acc);

%% Generate Measurements
disp('Generating measurements.') 

%Assumed standard deviations of the measurement noise components.
sigmaR = 10;
sigmaAz = 0.1*deg;
%Square root measurement covariance matrix; assume no correlation.
SR = diag([sigmaR,sigmaAz]);

[zMeasCart, SRMeasCart, tMeas, zMeasPol] = genMmts( ...
    ZCartTrue, ZPolTrue, PD, lambdaV, SR, mmt_space, Ts);
save("swarm_mmts", "zMeasCart", "SRMeasCart", "tMeas")

%% Motion Model

% Position space dimension
spaceDim = 2;

% State motion model
motion_models = [
    "Gauss-Markov Velocity (GMV)" 
    "Nearly Constant Velocity (NCV)"
    "Nearly Constant Acceleration (NCA)"];
% Extra parameters required by motion models
%   NCV: xDim = 4
%   - maximum acceleration: maxAcc (m/s/s)
%   NCA: xDim = 6
%   - maximum jerk: maxJrk (m/s/s/s)
%   GMV: xDim = 4
%   - maximum acceleration: maxAcc (m/s/s)
%   - maneuver decorrelation time: tau (s)
maxAcc = 9.8*2; % in g's
maxJrk = 0.1;
tau = 25*Ts;   % maneuver decorrelation time (s)
xDim = 4;
motion_model = Tracker.MotionModel( ...
    "Type", "GMV", "StateDim", xDim, ...
    "SpaceDim", spaceDim, "RevisitTime", Ts, ...
    "MaxKinVar", maxAcc, "Tau", tau);
% State motion model method - see displayTracksSwarm.m
motion_model_method = motion_model.Type;

%% State Initialization
init_methods = ["One-Point", "Two-Point", "Three-Point Heuristic"];
state_init = Tracker.StateInitialization(...
    "Type", "Two-Point", "VelMin", Vmin, "VelMax", Vmax, ...
    "StateDim", xDim, "SpaceDim", spaceDim, ...
    "MotionModelType", motion_model.Type);
% State initialization method - see displayTracksSwarm.m
init_method = state_init.Type;

%% Track Filter
%Now for the tracker with integrated track initiation/ termination.
disp('Running the integrated tracking algorithm.') 

%We will save the value of each track at each time. Thus, a structure
%holds the track states at each time. The first field in state_st are the
%actual values, the fourth field in state_st is an ID so that we can
%associate the tracks across time.

state_st = struct('x',[],'S',[],'r',[],'ID',[]);
state_st = repmat(state_st, [numSamples,1]);

track_st = struct('x',[],'S',[],'r',[],'ID',[],'scan_num',[],'num_hits',[]);

clear two_point_init
for curScan = 1:numSamples
    tCur = tMeas{curScan};
    zCurCart = zMeasCart{curScan};
    zPolCur = zMeasPol{curScan};
    SRCurCart = SRMeasCart{curScan};

    numMeas = size(zCurCart,2);

    % State initialization
    [xNew, SNew] = state_init.Initialize(tCur, zCurCart, SRCurCart);
    numStates = size(xNew, 2);

    fprintf("Scan No. %d\n", curScan)
    fprintf("Mmts: %d, States: %d\n", numMeas, numStates)

    %Initialization existence probabilities
    rNew = PIsRealInit*ones(1,numStates);
    
    %Generate a UUID for each potential track so that we can associate
    %tracks over time to draw lines for display. The UUIDs are 36-character
    %strings. For simplicity, we are just using the hash values of the
    %UUIDs so that they can be easily compared with >, = ,< for sorting.
    IDNew = zeros(1,numStates);
    for curNewTrack = 1:numStates
        [~,IDNew(curNewTrack)] = genUUID();
    end
    
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
                x(:,curTar),S(:,:,curTar), ...
                motion_model.F, ...
                motion_model.SQ);
        end
        %Update the target existence probabilities with the Markov 
        %switching model.
        r = PStayAlive*r;

        % Calculate mmt Jacobians
        measJacobDets = zeros(1,numMeas);
        for curMeas = 1:numMeas
            % measJacobDets(curTar) = det(calcPolarConvJacob( ...
            %     zPolCur(:,curTar),0,true));
            measJacobDets(curMeas) = det(calcPolarJacob( ...
                zCurCart(:,curMeas),0,true));
        end

        %The inclusion of r takes into account the track existence
        %probabilities.
        cm_obj.ClutterMap(zCurCart, curScan);
        lambda = cm_obj.Map(cm_obj.MmtIdx);
        [A,xHyp,PHyp] = makeStandardCartOnlyLRMatHyps(x,S,zCurCart,SRCurCart,[], ...
            PD,lambda,r,gammaVal,[]);
        
        [xUpdate,PUpdate,rUpdate,probNonTargetMeas] ...
            = singleScanUpdateWithExistence(xHyp,PHyp,PD,r,A, ...
            algSel1,algSel2,param3);
        
        %Determine which tracks to drop because their existence
        %probabilities are below the termination probability.
        sel = rUpdate>PTerminate;
        numTargetsCur = sum(sel);
        xUpdate = xUpdate(:,sel);
        PUpdate = PUpdate(:,:,sel);
        rUpdate = rUpdate(sel);
        xID = xID(sel);

        SUpdate = zeros(size(PUpdate));
        for curTar = 1:numTargetsCur
            SUpdate(:,:,curTar) = cholSemiDef(PUpdate(:,:,curTar),'lower');
        end
        
        if length(probNonTargetMeas) == length(rNew)
            rNew = probNonTargetMeas'.*rNew;
        end
        
        sel = rNew>PTerminate;
        xNew = xNew(:,sel);
        SNew = SNew(:,:,sel);
        rNew = rNew(sel);
        IDNew = IDNew(sel);
        
        state_st(curScan).x = [xNew, xUpdate];
        state_st(curScan).ID = [IDNew, xID];
        state_st(curScan).S = cat(3,SNew,SUpdate);
        state_st(curScan).r = [rNew, rUpdate];
    else
        state_st(curScan).x = xNew;
        state_st(curScan).ID = IDNew;
        state_st(curScan).S = SNew;
        state_st(curScan).r = rNew;
    end

    %% Populate track_st
    track_st = get_tracks(curScan, state_st(curScan), track_st);

end

%% Plots

% Estimated clutter map
cm_obj.PlotCLutterMap

% Tracks
displayTracksSwarm

%% LICENSE:
%
%The source code is in the public domain and not licensed or under
%copyright. The information and software may be used freely by the public.
%As required by 17 U.S.C. 403, third parties producing copyrighted works
%consisting predominantly of the material produced by U.S. government
%agencies must provide notice with such work(s) identifying the U.S.
%Government material incorporated and stating that such material is not
%subject to copyright protection.
%
%Derived works shall not identify themselves in a manner that implies an
%endorsement by or an affiliation with the Naval Research Laboratory.
%
%RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
%SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY THE NAVAL
%RESEARCH LABORATORY FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS
%OF RECIPIENT IN THE USE OF THE SOFTWARE.
