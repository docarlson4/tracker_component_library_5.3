function [ZCart, ZPol, Ts, num_tgt, max_acc] = genSwarm(Ts, Ns, Rmin, Rmax, Vmin, Vmax, lam)
%GENSWARM Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT
%   Ts - sample time or inverse rotation rate
%   Ns - nmber of samples or scans or revisits
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

if nargin == 0
%     rng(0)
    addpath c:\Users\dougo\OneDrive\Documents\MATLAB\
    Constants
    Ts = 12;
    Ns = 100;
    Rmin = 5*km;
    Rmax = 40*km;
    Vmin = 25*mph;
    Vmax = 600*mph;
    lam = 5;
    genSwarm(Ts, Ns, Rmin,Rmax, Vmin, Vmax, lam);
end

%%
f = @(x) fTransCoordTurn2D(Ts,x);

pois = PoissonD;

num_tgt = pois.rand(1, lam);
ZCart = cell(Ns,1);
ZPol = cell(Ns,1);

% Polar location of targets in NED
% rng_tgt = Rmin + (Rmax-Rmin)*rand(num_tgt, 1);
rng_tgt = Rmax*ones(num_tgt, 1);
alp_tgt = 2*pi*rand(num_tgt, 1);
pos = [rng_tgt.*sin(alp_tgt), rng_tgt.*cos(alp_tgt)]';

% Velocity of targets
speed = Vmin + (Vmax-Vmin)*rand(num_tgt,1);
bet_tgt = alp_tgt + (-pi/4 + pi/2*rand(num_tgt,1)) + pi;
vel = [speed.*sin(bet_tgt), speed.*cos(bet_tgt)]';

% Angular speed
rv = cross([pos;zeros(1,num_tgt)],[vel;zeros(1,num_tgt)],1);
w = rv(3,:)./(rng_tgt.*rng_tgt)';
% w = w .* (2 * (2*rand(1,num_tgt)-0));
rCtr = speed./w';
max_acc = rCtr.*w'.^2;

% Target lifetime
num_lt = 3 + pois.rand([num_tgt,1],Ns/3-3);
xCart = cell(num_tgt, 1);
xPol = cell(num_tgt, 1);
xDot = cell(num_tgt, 1);
x0 = [pos;vel;w];
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

% First detection
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
