function triplets = FindConnectedPaths(A21, A32)
%FINDCONNECTEDPATHS Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-06-07 10:30

%%
[n2, ~] = size(A21);
[~, n2b] = size(A32);

assert(n2 == n2b, 'Incompatible dimensions between A21 and A32');

triplets = [];

% Loop over all possible (i, j, k) triplets
for j = 1:n2
    from1 = find(A21(j, :));  % i's connected to j
    to3   = find(A32(:, j));  % k's connected from j

    for i = from1
        for k = to3'
            triplets = cat(1, triplets, [i,j,k]);
        end
    end
end

