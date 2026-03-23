function [xnodes, wnorm] = GaussHermiteNorm(n)
%GAUSSHERMITENORM  n-point Gauss-Hermite quadrature for N(0,1)
%
%   Approximates expectation under the standard normal distribution:
%       E[f(X)] = ∫ f(x) N(x;0,1) dx
%   ≈ sum_j wnorm(j) * f(xnodes(j))
%
% INPUT
%   n - number of quadrature points
%
% OUTPUT
%   xnodes - nodes for evaluation (size n×1)
%   wnorm  - corresponding normalized weights (size n×1)
%
% USAGE
%   [xnodes, wnorm] = GaussHermiteNorm(n)
%

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-09-20 16:07

if nargin == 0
    % Example: compute E[X^2] under N(0,1), which should be 1
    [xnodes, w] = GaussHermiteNorm(5);   % 5-point rule
    expectation = sum(w .* (xnodes.^2));
    disp(expectation)
    return
    % Output: 1.000000000000
end

%--- Construct Hermite polynomial coefficients H_n
% Physicists' Hermite polynomials via recurrence
H = SeaClutter.Hermite(n);

%--- Roots of H_n(y)
rootsHn = sort(roots(H{n+1}));

%--- Weights for exp(-y^2) quadrature
w = zeros(n,1);
for i = 1:n
    y = rootsHn(i);
    Hnm1 = polyval(H{n}, y); % H_{n-1}(y_i)
    w(i) = (2^(n-1) * factorial(n) * sqrt(pi)) / (n^2 * (Hnm1^2));
end

%--- Scale to standard normal
xnodes = sqrt(2) * rootsHn;   % nodes for N(0,1)
wnorm  = w / sqrt(pi);        % normalized weights

