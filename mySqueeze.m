function [x] = mySqueeze(x)
%MYSQUEEZE Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-07-20 14:56

%%
z = size(x);

x = reshape(x,[z(2:end),1]);


end
