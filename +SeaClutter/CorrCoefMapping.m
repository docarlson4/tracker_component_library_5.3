function [Rout, Rgaus] = CorrCoefMapping(params, num_gauss, MC_INTEG)
%CORRCOEFMAPPING Summary of this function goes here
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

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-09-20 17:22

%%

if nargin < 2
    num_gauss = 1001;
end
if nargin < 3
    MC_INTEG = false;
end

if MC_INTEG
    rng(1234); % for reproducibility
    ns = 10000;
    xs = randn(ns,1);
else
    ns = 20;
    [xn, wnorm] = SeaClutter.GaussHermiteNorm(ns);
end

alpha = params.alpha;
beta = params.beta;
switch params.type_out
    case 'Gamma'
        Q = @(x,alpha,beta) GammaD.invCDF(erfc(x)/2,alpha,beta);
    case 'Pareto'
        Q = @(x,alpha,beta) ParetoTypeID.invCDF(erfc(x)/2,alpha,beta);
    otherwise
end

Rgaus = linspace(-1,1,num_gauss)';
num_alp = numel(alpha);
Rout = zeros(num_gauss, num_alp);
for kN = 1:num_alp

    fs = ones(1,6);
    for n = 1:5
        if MC_INTEG
            arg = SeaClutter.Hermite(n+1, xs/sqrt(2)) .* Q(xs/sqrt(2),alpha(kN),beta(kN));
            fs(n+1) = mean(arg);
        else
            arg = SeaClutter.Hermite(n+1, xn/sqrt(2)) .* Q(xn/sqrt(2),alpha(kN),beta(kN));
            fs(n+1) = wnorm' * arg;
        end
        % figure
        % plot(xs,Hermite(n+1, xs/sqrt(2)),'.')
        % drawnow
    end
    n = 0:5;
    fs = (fs.^2)./(2.^n.*factorial(n));

    Rout(:,kN) = alpha(kN) * (Rgaus.^n * fs' - 1);
end




