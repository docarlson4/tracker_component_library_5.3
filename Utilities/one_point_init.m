function [xNew, SNew] = one_point_init(zCur, SRCur, Vmax, xDim)
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

[zDim, numMeas] = size(zCur);

%First, we create potential tracks for all of the observations.
xNew = zeros(xDim,numMeas);
SNew = zeros(xDim,xDim,numMeas);
xNew(1:zDim,:) = zCur;
SNew(1:zDim,1:zDim,:) = SRCur;
%The uncertainty for the unknown velocity
for curDim = (zDim+1):xDim
    SNew(curDim,curDim,:) = Vmax/sqrt(2+zDim);
end


end
