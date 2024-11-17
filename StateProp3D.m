function [xList] = StateProp3D(xInit, times, aDrift)
%STATEPROP3D Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT
%   xInit - initial state vector [posECF;velENU;angularSpeedUp]
%   times - sample times for evaluation, 1 x Ns
%   aDrift - drift function per reference [1]
%
% OUTPUT
%   xInit - final state vector [posECF;velENU;angularSpeedUp] x Ns
%
% USAGE
%   [xList] = StateProp3D(xInit, times, aDrift)
%
% REFERENCE
% [1] D. F. Crouse, "Simulating aerial targets in 3D accounting for the
%     Earth's curvature," Journal of Advances in Information Fusion, vol.
%     10, no. 1, Jun. 2015.

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-11-03 16:34

%The initial local coordinate system is defined as East-North-Up,
%corresponding to the x, y and z axes. This lets us easily translate the
%initial heading given by the indirectGeodeticProb function into components
%for the initial heading. The basis vectors at the starting point are
pos_tgt_ecf = xInit(1:3,1);
pos_tgt_geo = Cart2Ellipse(pos_tgt_ecf);
uInit = getENUAxes(pos_tgt_geo);

%When traveling in a straight line on an ellipsoidal Earth, the local basis
%vectors evolve according to
uDyn=@(u,x,t)uDotEllipsoid(u,x);

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
xList = RungeKCurvedAtTimes(xInit, uInit, times, aDrift, uDyn);

end
