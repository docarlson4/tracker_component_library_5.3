clc
clear
close all

cd(fileparts(which(mfilename)))
restoredefaultpath

dbstop if error

Constants

% Doug Carlson (doug.o.carlson@gmail.com), 2024-04-21 14:41

%% Parameters
% rng(1110)

%% Square
wid = 3;
hgt = 2;

xc = 1;
yc = 1;

% [x0,y0] = GenSquare(xc,yc,wid,hgt);
[x0,y0] = GenAnnularSector(xc,yc,0.5,22*deg);

%% Measurement
xm = randn;
ym = randn;

%% Create Inflated or Deflated Shape
rat = wid/hgt;
if abs(xm - xc) < abs(ym - yc) * rat
    hgt_m = abs(ym - yc) * 2;
    wid_m = hgt_m * rat;
else
    wid_m = abs(xm - xc) * 2;
    hgt_m = wid_m / rat;
end
[x_m,y_m] = GenSquare(xc,yc,wid_m,hgt_m);
vol_m = wid_m * hgt_m;

%% Plots
figure('Position',[950,200,560,420])
hold on, grid on, box on
plot(x0,y0,'-',DisplayName="Original",LineWidth=2)
plot(xm,ym,'o',DisplayName="Measurement",LineWidth=2)
plot(x_m,y_m,'--',DisplayName="New",LineWidth=2)
axis([-1 1 -1 1]*5)
axis equal
legend

%% Helper Functions
function [x,y] = GenAnnularSector(xc,yc,dr,da)
rc = sqrt(xc^2 + yc^2);
ac = atan2(yc,xc);
r1 = rc-dr/2; r2 = rc+dr/2;
a1 = ac-da/2; a2 = ac+da/2;
dela = da/11;
a = a1:dela:a2;
fa = fliplr(a);
x = [r2*cos(a), r1*cos(fa), r2*cos(a1)] + xc;
y = [r2*sin(a), r1*sin(fa), r2*sin(a1)] + yc;
end

function [x,y] = GenSquare(xc,yc,w,h)
x = [-1 -1 1 1 -1]*w/2 + xc;
y = [-1 1 1 -1 -1]*h/2 + yc;
end

