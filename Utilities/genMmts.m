function [zMeasCart, SRMeasCart, zMeasPol] = genMmts( ...
    ZCartTrue, ZPolTrue, PD, lambdaV, mmt_space, radar_obj, numStnry)
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

if nargin < 7
    numStnry = 0;
end
%Assumed standard deviations of the measurement noise components.
sigmaR = radar_obj.RangeUnc;
sigmaAz = radar_obj.AzimuthUnc;
sigmaRr = radar_obj.RangeRateUnc;
%Square root measurement covariance matrix; assume no correlation.
SR = diag([sigmaR,sigmaAz,sigmaRr]);

Ts = radar_obj.FrameTime;

AmbRR = radar_obj.RangeRateAmb;

%%
maxRange = mmt_space(2,1);

% Stationary clutter
km = 1e3;
rs = 0*km + 5*km * sqrt(rand(1,numStnry));
azs = rand(1,numStnry) * 2*pi;
rrs = 0.1*randn(1,numStnry);
% Compute the wrapped range-rate values
rrs = mod(rrs + AmbRR/2, AmbRR) - AmbRR/2;
ts = azs/(2*pi)*Ts;
IDs = 2*ones(1,numStnry);

%Cubature points for measurement conversion.
zDim = size(SR,1);
[xi,w] = fifthOrderCubPoints(zDim-1);

numSamples = length(ZCartTrue);

%Generate measurements and false alarms for each scan.
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
    zPolCur = zeros(zDim+2,numMeas);
    curDet = 1;

    if(numMeas>=0)
        %Generate the detection from the targets, if any.
        for curTar = 1:numTargets
            if(isDet(curTar)) && ( ZPolTrue{curScan}(1,curTar) <= maxRange)
                zPolCur(1:3,curDet) = ZPolTrue{curScan}(1:3,curTar) ...
                    + SR*randn(zDim,1);
                % Compute the wrapped range-rate values
                zPolCur(3,curDet) = mod(zPolCur(3,curDet) + AmbRR/2, AmbRR) - AmbRR/2;
                az = zPolCur(2,curDet);
                zPolCur(4,curDet) = ( az/(2*pi) + (curScan-1) ) * Ts;
                zPolCur(5,curDet) = 0;
                curDet = curDet+1;
            end
        end

        %Generate the false alarm detections, if any. 
        rClutBounds = mmt_space(:,1);%[minRange;maxRange];
        azBounds = mmt_space(:,2);%[-pi;pi];
        rrBounds = mmt_space(:,3);%
        for curFalse = 1:numFalse
            %r = UniformD.rand(1,rClutBounds);
            xDim = size(rClutBounds,2);
            r = (rClutBounds(2,:)-rClutBounds(1,:))' ...
                .* rand(xDim, 1) + rClutBounds(1,:)';
            az = UniformD.rand(1,azBounds);
            rr = UniformD.rand(1,rrBounds);
            % Compute the wrapped range-rate values
            rr = mod(rr + AmbRR/2, AmbRR) - AmbRR/2;

            t = ( az/(2*pi) + (curScan-1) ) * Ts;

            zPolCur(:,curDet) = [r; az; rr; t; 1];
            curDet = curDet+1;
        end
            
        %We will now convert the measurements into Cartesian coordinates as
        %we are using a converted-measurement filter.
        zPols = [rs;azs;rrs] + SR*randn(3,numStnry);
        zMeasPol{curScan} = [zPolCur, [zPols;ts + (curScan-1)*Ts;IDs]];
        [zMeasCart{curScan},RMeasCart] = pol2CartCubature( ...
            zMeasPol{curScan}(1:2,:),SR(1:2,1:2),0,true,[],[],[],xi,w);
        
        %Take the lower-triangular square root of the covariance matrices.
        for curMeas = 1:numMeas + numStnry
            RMeasCart(:,:,curMeas) = chol(RMeasCart(:,:,curMeas),'lower');
        end
        SRMeasCart{curScan} = RMeasCart;
    end
end



end
