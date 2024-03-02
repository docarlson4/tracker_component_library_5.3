function [state_st] = get_trajectory(traj_st, center)
%GET_TRAJECTORY Summary of this function goes here
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

% Developed in Matlab 9.11.0.1873467 (R2021b) Update 3 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2022-06-04 09:59

if nargin < 2
    center = false;
end

%% Constants
deg = pi/180;
km = 1e3;

%% Compass directions and angles in True North CS
compass_dir = {'N','NNE','NE','ENE','E','ESE','SE','SSE','S','SSW','SW','WSW','W','WNW','NW','NNW'};
compass_ang = (0:22.5:359);

%%
T = traj_st.T;
v0 = traj_st.v0;

if ~isfield(traj_st, 'pos0')
    traj_st.pos0 = [0;0;0];
end

num_leg = numel(traj_st.traj_dir);
num_seg = 2*num_leg-1;

modes = false(1,num_seg); modes(2:2:num_seg) = true;
times = nan(1,num_seg); times(1:2:num_seg) = traj_st.traj_time;
%acc = zeros(1,num_seg); acc(2:2:num_seg) = traj_st.traj_acc;

ang = zeros(1,num_leg);
for kL = 1:num_leg
    ang(kL) = compass_ang(matches(compass_dir, traj_st.traj_dir(kL)));
end
del_phi = -diff(ang)*deg;

for kA = 1:num_leg-1

    if traj_st.traj_acc(kA) > 0 && del_phi(kA) < 0
        del_phi(kA) = del_phi(kA) + 2*pi;
    elseif traj_st.traj_acc(kA) < 0 && del_phi(kA) > 0
        del_phi(kA) = 2*pi - del_phi(kA);
    end
end

N_turn = ceil( abs(del_phi)*v0./(abs(traj_st.traj_acc)*T) );
N = times/T; N(2:2:num_seg) = N_turn;
%times(2:2:num_seg) = N_turn*T;
w = zeros(1,num_seg); w(2:2:num_seg) = sign(traj_st.traj_acc).*abs(del_phi)./(T*N_turn);
N = ceil(N);
num_scan = sum(N)+1;
N_cum = [0, cumsum(N)]+1;
pos = traj_st.pos0;
heading = [sin(ang(1)*deg);cos(ang(1)*deg);0];
vel = v0*heading;

s = [pos;vel;0]*ones(1,num_scan);
t = zeros(1,num_scan);
mode = zeros(1,num_scan);
idx = [1:2, 4:5, 7];
for kS = 2:num_seg+1

    for k = N_cum(kS-1)+1:N_cum(kS)
        s(7,k-1) = w(kS-1);
        s(idx,k) = FCoordTurn2D(T,s(idx,k-1)) * s(idx,k-1);
        t(k) = (k-1)*T;
        %mode(k) = modes(kS-1);
        mode(k) = sign(w(kS-1));
    end
end
if center
    s(1:2,:) = s(1:2,:) - (max(s(1:2,:),[],2) + min(s(1:2,:),[],2))/2;
end

state_st.s = s;
state_st.mode = mode;
state_st.t = t;


end
