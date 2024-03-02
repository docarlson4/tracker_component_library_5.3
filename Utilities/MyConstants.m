% Physical and mathematical constants

% angles
d2r = pi/180;
r2d = 180/pi;
deg = pi/180;
rad = 180/pi;
mrad = 1e-3;

% Frequency
Hz = 1;
THz = 1e12;
GHz = 1e9;
MHz = 1e6;
kHz = 1e3;

% Impedance
MOhm = 1e6;
kOhm = 1e3;
Ohm = 1;
S = 1/Ohm; % Siemens

% E-M
c0 = 299792458;
mu0 = 4e-7*pi;
ep0 = 1/(mu0*c0^2);
et0 = sqrt(mu0/ep0);
eV = 1.602176487e-19; % J
keV = 1e3*eV;
MeV = 1e6*eV;

% Electrical
pF = 1e-12;   % pico-Farad
uF = 1e-6;    % micro-Farad
mH = 1e-3;    % milli-Henry
uH = 1e-6;
qE = 1.60217662e-19; % Coulomb

% Power
W = 1;
Watts = 1;
kW = 1e3;
MW = 1e6;

% Thermodynamics
kB = 1.38064852e-23; % Joules/Kelvin

% Time
second = 1;
ms = 1e-3;
us = 1e-6;
ns = 1e-9;

% Distance
km = 1e3;
cm = 1e-02;
mm = 1e-03;
um = 1e-06;
nm = 1e-09;
pm = 1e-12;
fm = 1e-15;
inches = 0.0254;
mils = inches*1e-3;
ft = 12*inches;
miles = 5280*ft;

% Mass
kg = 1;
gm = 1e-3;

% Quantum
hbarc = 197.327053*MeV*fm;
hbar = 1.05457266e-34; % J-s
me = 510998.909848 * eV / c0^2; % kg
Angstrom = 1e-10;

% Speed
mps = 1;
kps = 1e3*mps;
mph = 5280*12*0.0254/3600;
knots = 0.51444444; % 463/900
kph = 1/3.6;

% Earth
RE = 6378137;
RP = 6356752.3142;
REmean = 6371008.7714;
ME = 5.972e24;
g0 = 9.80665;
omegaEarth = 7.2921158553e-5;

% Mars
RM = 3396.19 * km;
muM = 42828.3 * km^3/second^2;
gM0 = 9.80665;
omegaMars = 7.09003e-5;

% Useful units
charAngstrom = char(197);
RPM = 2*pi/60;

% Numerical Constants
eulerGamma = exp(0.577215664901532860606512090082402431042159335);
gammaE = eulerGamma;

% Gravity
Ggrav = 6.67430e-11;





