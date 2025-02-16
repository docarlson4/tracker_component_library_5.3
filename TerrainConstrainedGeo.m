function [zTgt] = TerrainConstrainedGeo(z, hp, zRx, AC, maxIter, tol)
%TERRAINCONSTRAINEDGEO geolocates a target based on range-direction cosine
%mmts in the antenna frame of an airborne platform, using Newton's method
%
% INPUT
%   z - [r;u]: 2 x N range-direction cosine mmts
%   hp - height of the platform
%   zRx - position of platform in ECEF
%   AC - transform from ECEF to antenn coordinates ACS
%   maxIter - limit the number of iterations of Newton's method
%   tol - limit update tolerance of Newton's method
%
% OUTPUT
%   zTgt - target location in ECEF
%
% USAGE
%   [zTgt] = TerrainConstrainedGeo(z, hp, zRx, AC, maxIter, tol)
%

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-02-15 13:48

%% Loop Control
if nargin < 5
    maxIter = 10;
end

if nargin < 6
    tol = 1e-9;
end

%% Constants
a = Constants.WGS84SemiMajorAxis;
f = Constants.WGS84Flattening;
b = a*(1-f);
aE = (2*a+b)/3;
aabInv = diag(1./[a,a,b]);

%% Initial guess
r = z(1,:);
u = z(2,:);
numC = size(u,2);

v0 = (aE*aE - r.*r - (aE+hp).*(aE+hp))./(2*(aE+hp).*r);
w0 = sqrt(1 - u.*u - v0.*v0);

uvw = [u;v0;w0];
zTgt = zRx + r .* (AC' * uvw);

%% Newton method loop
iter = 0;
delv = inf;
while (iter < maxIter) && (max(abs(delv)) > tol)
    iter = iter+1;
    
    yTgt = aabInv*zTgt;
    f = sum(yTgt.*yTgt,1) - 1;
    dzTgt = r.* (AC' * [zeros(1,numC); ones(1,numC); -v0./w0]);
    dyTgt = aabInv*dzTgt;
    df = 2*sum(yTgt.*dyTgt, 1);

    delv = f./df;
    v0 = v0 - delv;
    w0 = sqrt(1 - u.*u - v0.*v0);

    uvw = [u;v0;w0];
    zTgt = zRx + r .* (AC' * uvw);

end

% disp(iter)
% disp(max(abs(delv)))

