clc
clearvars
close all

cd(fileparts(which(mfilename)))

startup_tcl

dbstop if error

MyConstants

% Doug Carlson (doug.o.carlson@gmail.com), 2024-04-10 19:20

rng(0)

%% Monte Carlo Samples
numMC = 1;

%% Measurement Model
zDim = 2;
mmt_space = [
    0, 0
    1e3, 400];
mmt_vol = prod(diff(mmt_space)); % m^2

lambda = 1e-4; % 1/scan/m^2
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
PG = 0.9999;
gammaVal = ChiSquareD.invCDF(PG,2);


%% Truth Data
xDim = 4;
s0 = [130;200;35;0] * ones(1,numSamples);
for curScan = 2:numSamples
    s0(:,curScan) = F*s0(:,curScan-1);
end

%% Measurements
zMeasCart = cell(numSamples,numMC);
for curMC = 1:numMC

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
            zMeasCart{curScan, curMC} = zCartCur;
        end

    end

end

%% Filter
state_st = struct('x',[],'SP',[],'r',[],'ID',[]);
state_st = repmat(state_st, [numSamples,numMC]);
   
for curMC = 1:numMC
%We will save the value of each track at each time. Thus, a structure
%holds the track states at each time. The first field in state_st are the
%actual values, the fourth field in state_st is an ID so that we can
%associate the tracks across time.

clear two_point_init
x = zeros(4,0);
for curScan = 1:numSamples
    tCur = (curScan-1)*T;
    zCur = zMeasCart{curScan};
    SRCur = SR;
    numMeas = size(zCur,2);

    if isempty(x)
        % case "Two-Point" initialization
        [x, SP] = two_point_init(zCur(:,1), SRCur, tCur, [0,50], xDim);

        %% Initialize existence probabilities
        r = PIsRealInit;

        state_st(curScan, curMC).x = x;
        state_st(curScan, curMC).SP = SP;
        state_st(curScan, curMC).r = r;        
    else
        %% Predict the tracks to the current time.
        [xPred,SPPred] = sqrtDiscKalPred(x,SP,F,SQ);
        PPred = SPPred*SPPred';

        % Update the target existence probabilities with the Markov 
        % switching model.
        rPred = p11*r;

        %% Predicted mmt and validated innovation
        zPred = H*xPred;
        nu = zMeasCart{curScan} - zPred;
        S = H*PPred*H' + R;
        lgcl_gam = dot(nu,S\nu) < gammaVal;
        sprintf("Num in gate: %d", sum(lgcl_gam))
%         disp([dB(gammaVal), dB(dot(nu,S\nu))])
        nu = nu(:,lgcl_gam);

        %% Kalman Gain
        K = PPred * (H' * (S\eye(zDim)));

        %% Likelihoods: Lam
        gauss = GaussianD;
        Lam = gauss.PDF(nu,[],S);
        del = PD*PG * (1 - mmt_vol/lambdaV * sum(Lam,2));
        bet0 = (1-PD*PG)*rPred/(1 - del*rPred);
        beti = PD*PG * mmt_vol/lambdaV * rPred * Lam/(1 - del*rPred);
        r = sum([bet0,beti]);

        %% State Updates

        nub = nu * beti';
        x = xPred + K * nub;

        Pc = PPred - K*S*K';

        nunu = nub*nub';
        for k = 1:sum(lgcl_gam)
            nunu = nunu + beti(k)*(nu(:,k)*nu(:,k)');
        end
        Pw = K*(nunu)*K';

        P = bet0*PPred + (1-bet0)*Pc + Pw;
        SP = cholSemiDef(P,'lower');

        state_st(curScan, curMC).x = x;
        state_st(curScan, curMC).SP = SP;
        state_st(curScan, curMC).r = r;        
    end    

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

s = [state_st.x];
r = [state_st.r]
plot(s(1,:),s(2,:),'r*')

axis equal
% axis(mmt_space(:)')











