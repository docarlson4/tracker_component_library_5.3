function [h, H] = Hermite(n, x)
%HERMITE recursively generates Hermite polynomial coefficients if only "n"
%is supplied. If domain "x" is also supplied, then the polynomial for the
%largest "n" is provided.
%
%Note: the algorithm aborts at n = 263, as beyond this value nans are
%produced.
%
% INPUT
%   n - order of Hermite polynomial
%   x - domain
%
% OUTPUT
%   h - coefficients or Hn(x)
%   H - coefficients
%
% USAGE
%   [h, H] = Hermite(n, x)

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-09-20 15:40

%%

if nargin == 0
    n = 5;
    [h] = Hermite(n);
    disp(h)

    x = -2:0.01:2;
    figure
    hold on, grid on, box on
    for k = 1:n
        plot(x,Hermite(k,x), DisplayName="H_"+num2str(k-1)+"(x)", LineWidth=2)
    end
    ylim([-30,30])
    legend
end

if nargin == 1
    x = nan;
end

% Construct Hermite polynomial coefficients H_n Physicists' Hermite
% polynomials via recurrence
H = cell(n+1,1);
H{1} = 1;          % H_0(y) = 1
H{2} = [2 0];      % H_1(y) = 2y
for k = 2:n
    % H_k(y) = 2y H_{k-1}(y) - 2(k-1) H_{k-2}(y)
    H{k+1} = conv([2 0], H{k}) - 2*(k-1)*[0 0 H{k-1}];
    if any(isnan(H{k+1}))
        error("Reached numerical limit: " + num2str(k))
    end
end

h = H;

if nargin == 2
    h = polyval(H{n}, x);
end

if nargout == 0
    clear h
end
