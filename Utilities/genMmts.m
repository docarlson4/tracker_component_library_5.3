function [zMeasCart, SRMeasCart,tMeas, zMeasPol] = genMmts( ...
    ZCartTrue, ZPolTrue, PD, lambdaV, SR, mmt_space, Ts)
%GENMMTS Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT
%
%
% OUTPUT
%
%
% USAGE
%
%

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-23 18:44

%%
maxRange = mmt_space(2,1);

% Stationary clutter
km = 1e3;
numStnry = 20;
rs = 0*km + 5*km * sqrt(rand(1,numStnry));
azs = rand(1,numStnry) * 2*pi;

%Cubature points for measurement conversion.
zDim = size(SR,1);
[xi,w] = fifthOrderCubPoints(zDim);

numSamples = length(ZCartTrue);

%Generate measurements and false alarms for each scan.
tMeas = cell(numSamples,1);
zMeasCart = cell(numSamples,1);
zMeasPol = cell(numSamples,1);
SRMeasCart = cell(numSamples,1);
for curScan = 1:numSamples

    %Determine the number of false alarms to generate.
    numFalse = PoissonD.rand(1,lambdaV);
    numTargets = size(ZCartTrue{curScan},2);
    if(curScan <= numSamples)
        %Determine which, if any, targets should be detected.
        isDet = rand(numTargets,1) < PD;
    else
        %We are after the end; the tracks are now gone.
        isDet = [0;0];
    end

    %Allocate space for the detections.
    numMeas = numFalse+sum(isDet);
    zPolCur = zeros(2,numMeas);
    curDet = 1;

    if(numMeas>=0)
        %Generate the detection from the targets, if any.
        for curTar = 1:numTargets
            if(isDet(curTar)) && ( ZPolTrue{curScan}(1,curTar) <= maxRange)
                zPolCur(:,curDet) = ZPolTrue{curScan}(:,curTar) ...
                    + SR*randn(2,1);
                curDet = curDet+1;
            end
        end

        %Generate the false alarm detections, if any. 
        rClutBounds = mmt_space(:,1);%[minRange;maxRange];
        azBounds = mmt_space(:,2);%[-pi;pi];
        for curFalse = 1:numFalse
            %r = UniformD.rand(1,rClutBounds);
            xDim = size(rClutBounds,2);
            r = (rClutBounds(2,:)-rClutBounds(1,:))' ...
                .* rand(xDim, 1) + rClutBounds(1,:)';
            az = UniformD.rand(1,azBounds);

            zPolCur(:,curDet) = [r; az];
            curDet = curDet+1;
        end
            
        %We will now convert the measurements into Cartesian coordinates as
        %we are using a converted-measurement filter.
        zMeasPol{curScan} = [zPolCur, [rs;azs]];
        [zMeasCart{curScan},RMeasCart] = pol2CartCubature( ...
            zMeasPol{curScan},SR,0,true,[],[],[],xi,w);
        
        %Take the lower-triangular square root of the covariance matrices.
        for curMeas = 1:numMeas + numStnry
            RMeasCart(:,:,curMeas) = chol(RMeasCart(:,:,curMeas),'lower');
        end
        SRMeasCart{curScan} = RMeasCart;
        tMeas{curScan} = repmat((curScan-1)*Ts, 1, numMeas + numStnry);
    end
end



end
