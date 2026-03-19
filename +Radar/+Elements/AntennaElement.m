classdef (Abstract) AntennaElement < handle
%ANTENNAELEMENT  Abstract base class for phased-array element patterns.
%
%   All element models derive from this class and must implement
%   ElementGain(az, el).  Shared utilities (PeakGain, PlotPattern) are
%   provided here so every subclass gets them for free.
%
%   SUBCLASSES
%     CosineElement       – separable cos^(2N) model
%     PatchElement        – rectangular microstrip patch approximation
%     VivaldiElement      – tapered-slot (Vivaldi) approximation
%     CrossedDipoleElement– two orthogonal half-wave dipoles in quadrature
%     NumericalElement    – data-driven pattern via interpolation
%
%   All angles are in radians.

    properties
        Description (1,1) string = ""
    end

    % -----------------------------------------------------------------
    %  Abstract interface – every subclass must supply this
    % -----------------------------------------------------------------
    methods (Abstract)

        G = ElementGain(obj, az, el)
        %ELEMENTGAIN  Linear power gain at azimuth AZ, elevation EL.
        %
        %   AZ and EL must be conformable arrays (radians).
        %   G is a non-negative real array of the same size.
        %   Pattern is defined with broadside at (az=0, el=0) and the
        %   propagation vector convention:
        %     k_hat = [cos(el)*sin(az); sin(el); cos(el)*cos(az)]

    end

    % -----------------------------------------------------------------
    %  Shared concrete utilities
    % -----------------------------------------------------------------
    methods

        function [peakGain, peakAz, peakEl] = PeakGain(obj, azRange, elRange)
        %PEAKGAIN  Maximum element gain and its location on an az/el grid.
        %
        %   [G, az, el] = PeakGain()                 – 1° grid over ±89°
        %   [G, az, el] = PeakGain(azRange, elRange) – custom grids (rad)
            arguments
                obj
                azRange (1,:) double = linspace(-pi/2*0.989, pi/2*0.989, 179)
                elRange (1,:) double = linspace(-pi/2*0.989, pi/2*0.989, 179)
            end
            [az, el] = meshgrid(azRange, elRange);
            G        = obj.ElementGain(az, el);
            [peakGain, idx] = max(G(:));
            [iEl, iAz]      = ind2sub(size(G), idx);
            peakAz = azRange(iAz);
            peakEl = elRange(iEl);
        end

        function ax = PlotPattern(obj, az0, el0)
        %PLOTPATTERN  Surface plot of the element pattern in az/el space.
        %
        %   ax = PlotPattern()
        %   ax = PlotPattern(az0, el0)   – degree vectors
            arguments
                obj
                az0 (1,:) double = -89:89
                el0 (1,:) double = -89:89
            end
            deg      = pi/180;
            [az, el] = meshgrid(az0*deg, el0*deg);
            G        = obj.ElementGain(az, el);
            gainDB   = 10*log10(max(G, 10^(-50/10)));

            figure;
            surf(az0, el0, gainDB);
            shading flat;
            xlabel('\bf Azimuth (deg)');
            ylabel('\bf Elevation (deg)');
            peakDB = max(gainDB(:));
            title({ sprintf('\\bf\\fontsize{14} %s', obj.Description), ...
                    sprintf('\\bf\\fontsize{12} Peak = %.1f dB', peakDB) });
            ax = gca;
        end

    end
end
