function [xNew, SNew] = OnePoint(obj, zCur, SRCur)
%ONE_POINT_INIT Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-30 15:01

%%
Vmax = obj.VelMax;
Amax = obj.AccMax;
xDim = obj.StateDim;
spaceDim = obj.SpaceDim;
[zDim, numMeas] = size(zCur);
idxPos = 1:spaceDim;
idxVel = (idxPos) + spaceDim;
idxAcc = (idxVel) + spaceDim;

%First, we create potential tracks for all of the observations.
xNew = zeros(xDim,numMeas);
SNew = zeros(xDim,xDim,numMeas);
xNew(1:zDim,:) = zCur;
SNew(1:zDim,1:zDim,:) = SRCur;
%The uncertainty for the unknown velocity
for curDim = idxVel
    SNew(curDim,curDim,:) = Vmax/sqrt(2+zDim);
end
if obj.MotionModelType == "NCA"
    for curDim = idxAcc
        SNew(curDim,curDim,:) = Amax/sqrt(2+zDim);
    end
end


