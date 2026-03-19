clc
clear
close all

cd(fileparts(which(mfilename)))

dbstop if error

Constants

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-07-03 20:24

%%
x0 = 0; x1 = 1;
y0 = 0; y1 = 1;
xb = 0; yb = 1;
xm = (x0 + x1)/2; ym = (y0 + y1)/2;
m = (y1 - y0)/(x1 - x0);

dist = abs(m.*(xb-xm) - (yb-ym))./sqrt(m.^2+1)
xi = (xb + m.*(m.*xm+yb-ym))./(m.^2+1)
yi = (ym + m.*(m.*yb+xb-xm))./(m.^2+1)

lgcl = sqrt((xi - x0)^2 + (yi - y0)^2) <= sqrt((x0 - x1)^2 + (y0 - y1)^2) & ...
    sqrt((xi - x1)^2 + (yi - y1)^2) <= sqrt((x0 - x1)^2 + (y0 - y1)^2)
