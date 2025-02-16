function [zTgtECF, SPTgtECF, zTgtGEO, SPTgtGEO] = ...
    TerrainConstrainedGeoCubature(z, SR, hp, zRx, AC, maxIter, tol)
%TERRAINCONSTRAINEDGEOCUBATURE Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-02-16 08:34

%% Loop Control
if nargin < 6
    maxIter = 10;
end

if nargin < 7
    tol = 1e-9;
end

%% Update Input Sizes
[dimZ, numMmt] = size(z);
if size(SR,3) == 1
    SR = repmat(SR, 1,1,numMmt);
end
if size(hp,2) == 1
    hp = repmat(hp, 1,numMmt);
end

zTgtECF = zeros(3,numMmt);
SPTgtECF = zeros(3,3,numMmt);
zTgtGEO = zeros(3,numMmt);
SPTgtGEO = zeros(3,3,numMmt);

%% Cubature
[xi,w] = fifthOrderCubPoints(dimZ);
numCub = length(w);
for kM = 1:numMmt
    zCub = z(:,kM) + SR(:,:,kM) * xi;

    [zTgtEcfCub] = TerrainConstrainedGeo( ...
        zCub, hp(:,kM), zRx(:,kM), AC(:,:,kM), maxIter, tol);
    zTgtGeoCub = Cart2Ellipse(zTgtEcfCub);

    zTgtECF(:,kM) = zTgtEcfCub*w;
    zTgtGEO(:,kM) = zTgtGeoCub*w;
    for kC = 1:numCub
        dzCub = zTgtEcfCub(:,kC) - zTgtECF(:,kM);
        SPTgtECF(:,:,kM) = SPTgtECF(:,:,kM) + w(kC)*(dzCub*dzCub');
        dzCub = zTgtGeoCub(:,kC) - zTgtGEO(:,kM);
        SPTgtGEO(:,:,kM) = SPTgtGEO(:,:,kM) + w(kC)*(dzCub*dzCub');
    end
    SPTgtECF(:,:,kM) = cholSemiDef(SPTgtECF(:,:,kM),'lower');
    SPTgtGEO(:,:,kM) = cholSemiDef(SPTgtGEO(:,:,kM),'lower');
end

end
