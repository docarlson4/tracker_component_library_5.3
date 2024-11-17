function obj = genSwarm(obj)
% function [ZCart, ZPol, Ts, num_tgt, max_acc] = ...
%         genSwarm(Ts, Ns, Rmin, Rmax, Vmin, Vmax, lam)
%GENSWARM Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT
%   Ts - sample time or inverse rotation rate
%   Ns - number of samples or scans or revisits
%   Rmin, Rmax - min and max range
%   Vmin, Vmax - min and max speed
%   lam - average number of targets
%
% OUTPUT
%   x
%
% USAGE
%
%

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-02 13:36

%%

% Propagator
Ns = obj.ScanNumber;
Ts = obj.RadarObject.FrameTime;
times = (0:Ns-1)*Ts;

pois = PoissonD;

num_tgt = pois.rand(1, obj.TargetNumber);

ZCart = cell(Ns,1);
ZPol = cell(Ns,1);

% Polar location of targets, azimuth in NED
rng_tgt = diff(obj.RangeLimits)*ones(num_tgt, 1) ...
    + obj.RangeLimits(1);
eps_tgt = 15*obj.deg + 15*obj.deg * rand(num_tgt, 1);
alp_tgt = 2*pi*rand(num_tgt, 1);
pos_tgt_enu = [ ...
    rng_tgt.*sin(alp_tgt).*cos(eps_tgt), ...
    rng_tgt.*cos(alp_tgt).*cos(eps_tgt), ...
    rng_tgt.*sin(eps_tgt)]';
pos_tgt_ecf = obj.CE * pos_tgt_enu + obj.pos_ref_ecf;
pos_tgt_geo = Cart2Ellipse(pos_tgt_ecf);

% Velocity of targets
speed = diff(obj.VelocityLimits)*rand(num_tgt,1) ...
    + obj.VelocityLimits(1);
% Heading
bet_tgt = alp_tgt + (-pi/4 + pi/2*rand(num_tgt,1)) + pi;
vel_tgt_enu = [ ...
    speed.*sin(bet_tgt), ...
    speed.*cos(bet_tgt), ...
    zeros(num_tgt,1)]';

% Angular speed
rv = cross(pos_tgt_enu, vel_tgt_enu,1);
w = rv(3,:)./(rng_tgt.*rng_tgt)';
% w = w .* (2 * (2*rand(1,num_tgt)-0));
% rCtr = speed./w';
% max_acc = rCtr.*w'.^2;

%% Propagate 3D
%The continuous-time drift function for a 3D, flat-Earth dynamic model is
aFlatEarth3D=@(x,t)aCircTurn(x);

%The initial target state with GLOBAL position and LOCAL velocity is thus
xInit = [pos_tgt_ecf; vel_tgt_enu; w];

%The initial local coordinate system is defined as East-North-Up,
%corresponding to the x, y and z axes. This lets us easily translate the
%initial heading given by the indirectGeodeticProb function into components
%for the initial heading. The basis vectors at the starting point are
% uInit=getENUAxes(pos_tgt_geo);

%When traveling in a straight line on an ellipsoidal Earth, the local basis
%vectors evolve according to
uDyn=@(u,x,t)uDotEllipsoid(u,x);

%Perform Runge-Kutta integration to determine the location of the target.
%The function RungeKCurvedAtTimes maps the flat-Earth model to a curved
%Earth.
xList = zeros(7,Ns,num_tgt);
for kT = 1:num_tgt
    uInit = getENUAxes(pos_tgt_geo(:,kT));
    xList(:,:,kT) = RungeKCurvedAtTimes(xInit(:,kT), uInit, times, ...
        aFlatEarth3D, uDyn);
end
return
%% Target lifetime
num_lt = 3 + pois.rand([num_tgt,1],Ns/3-3);
xCart = cell(num_tgt, 1);
xPol = cell(num_tgt, 1);
xDot = cell(num_tgt, 1);
x0 = [pos_tgt_enu;vel_tgt_enu;w];
f = @(x) fTransCoordTurn2D(Ts,x);
for kT = 1:num_tgt
    xCart{kT} = x0(:,kT);
    xPol{kT} = [
        vecnorm( xCart{kT}(1:2,1) )
        atan2( xCart{kT}(1,1), xCart{kT}(2,1) )
        getRangeRate(xCart{kT}(1:4,1),true)];
    for kLT = 2:num_lt(kT)
        xCart{kT}(:,kLT) = f( xCart{kT}(:,kLT-1) );
        xDot{kT}(:,kLT) = (xCart{kT}(:,kLT) - xCart{kT}(:,kLT-1))/Ts;

        xPol{kT} = [xPol{kT}, [
            vecnorm( xCart{kT}(1:2,kLT) )
            atan2( xCart{kT}(1,kLT), xCart{kT}(2,kLT) )
            getRangeRate(xCart{kT}(1:4,kLT),true)] ];
    end
end

%% First detection
num_fd = zeros(num_tgt, 1);
for kFD = 1:num_tgt
    num_fd(kFD) = randi([1,Ns-num_lt(kFD)],1);
end

for kS = 1:Ns
    for kT = 1:num_tgt
        if num_fd(kT) == kS
            ZCart{kS} = xCart{kT}(:,1);
            ZPol{kS} = xPol{kT}(:,1);
        end
        if (num_fd(kT) < kS) && (kS <= num_fd(kT) + num_lt(kT))
            idx = kS - num_fd(kT);
            if xPol{kT}(1,idx) > Rmax
                continue
            end
            ZCart{kS} = [ZCart{kS}, xCart{kT}(:,idx)];
            ZPol{kS} = [ZPol{kS}, xPol{kT}(:,idx)];
        end
    end
end

end
