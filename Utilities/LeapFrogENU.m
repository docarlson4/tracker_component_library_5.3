function [sFinalENU, t] = LeapFrogENU(sInitENU, Ns, Ts, fProp)
%LEAPFROGENU Summary of this function goes here
%   Detailed explanation goes here
%
% INPUT
%   sInitENU - initial state vector [posENU; velENU]
%   Ns - number of samples
%   Ts - frame time
%   fProp - propagator function handle
%
% OUTPUT
%   sFinalENU - final state vector [posENU; velENU] x Ns
%   t - time samples, 1 x Ns
%
% USAGE
%   % Propagator
%   F = @(t) kron([1,t;0,1],eye(3));
%   fProp = @(t,s) F(t)*s;
%
%   [sFinalENU, t] = LeapFrogENU(sInitENU, Ns, Ts, fProp)
%

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-11-03 16:55

%%
sFinalENU = sInitENU * ones(1, Ns);
t = zeros(1, Ns);

% Azimuth
az = @(s) mod(atan2d(s(1,:),s(2,:)), 360);

az0 = 0;
for k = 1:Ns
    iter = 0;
    while iter < 5
        az1 = az(sFinalENU(:,k));
        dt = Ts*(az1-az0)/360;
        sFinalENU(:,k) = fProp(dt, sFinalENU(:,k));
        t(k) = t(k) + dt;
        az0 = az1;
        iter = iter+1;
    end
    if k < Ns
        sFinalENU(:,k+1) = sFinalENU(:,k);
        t(k+1) = t(k);
        az0 = az0-360;
    end
end

end
