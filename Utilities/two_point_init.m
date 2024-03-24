function [s, S] = two_point_init(zBuf, SBuf, tBuf, vLims)
%TWO_POINT_INIT Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-23 20:45

%%
s = [];
S = [];

%%
t1 = tBuf{1}(1,:);
x1 = zBuf{1}(1,:);
y1 = zBuf{1}(2,:);

t2 = tBuf{2}(1,:);
x2 = zBuf{2}(1,:);
y2 = zBuf{2}(2,:);

del_t = t2' - t1;
del_x = x2' - x1;
del_y = y2' - y1;
vx = del_x./del_t;
vy = del_y./del_t;
v = sqrt(vx.*vx + vy.*vy);

lgcl_v = (vLims(1) < v) & (v < vLims(2));

[i2,i1] = find(lgcl_v)
sx = x2(i2);
sy = y2(i2);
svx = zeros(size(sx));
svy = zeros(size(sy));
for k = 1:length(i1)
    svx(k) = vx(i2(k),i1(k));
    svy(k) = vy(i2(k),i1(k));
end
s = [sx;sy;svx;svy];

end
