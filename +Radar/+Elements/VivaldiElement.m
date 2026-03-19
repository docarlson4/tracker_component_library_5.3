classdef VivaldiElement < Radar.Elements.AntennaElement
%VIVALDIELEMENT  Approximate element pattern for a Vivaldi tapered-slot.
%
%   A Vivaldi element is a broadband, wide-scan aperture antenna whose
%   pattern is wider and flatter near broadside than a cosine^2 model.
%   The model used here is a half-power raised-cosine product:
%
%     G(az,el) = G0 · cos^(2·m)(az/2) · cos^(2·n)(el/2)
%
%   The half-angle argument produces a broader main lobe (≈ twice the
%   3-dB width of a cos^(2k) pattern) and a smoother roll-off, consistent
%   with measured Vivaldi characteristics at UHF/X-band.
%
%   CONSTRUCTION
%     e = VivaldiElement()
%     e = VivaldiElement(TaperExpHorz=1.5, TaperExpVert=2.0)
%
%   PROPERTIES
%     TaperExpHorz  m — horizontal taper exponent (default 1.5)
%     TaperExpVert  n — vertical taper exponent   (default 2.0)
%
%   Reference: D. Schaubert et al., "Wide-band Vivaldi arrays," 2003.
%   (Pattern approximation only; substrate and flare geometry not modelled.)

    properties
        TaperExpHorz (1,1) double {mustBePositive} = 1.5  % horizontal exponent m
        TaperExpVert (1,1) double {mustBePositive} = 2.0  % vertical exponent   n
    end

    properties (Access = private)
        NormFactor (1,1) double = 1
    end

    methods

        function obj = VivaldiElement(options)
            arguments
                options.TaperExpHorz (1,1) double {mustBePositive} = 1.5
                options.TaperExpVert (1,1) double {mustBePositive} = 2.0
            end
            obj.TaperExpHorz = options.TaperExpHorz;
            obj.TaperExpVert = options.TaperExpVert;
            obj.NormFactor   = obj.computeNorm();
            obj.Description  = sprintf( ...
                "Vivaldi element (m=%.4g, n=%.4g)", ...
                obj.TaperExpHorz, obj.TaperExpVert);
        end

        function G = ElementGain(obj, az, el)
        %ELEMENTGAIN  Half-angle raised-cosine pattern, normalised to 4π sr.
        %
        %   G is zero outside the front hemisphere (|az|, |el| > π/2).
            m = obj.TaperExpHorz;
            n = obj.TaperExpVert;
            G = obj.NormFactor .* cos(az/2).^(2*m) .* cos(el/2).^(2*n);
            G(abs(az) > pi/2 | abs(el) > pi/2) = 0;
        end

    end

    % -----------------------------------------------------------------
    methods (Access = private)
    % -----------------------------------------------------------------

        function nf = computeNorm(obj)
        %COMPUTENORM  Scale so total radiated power = 4π steradians.
            N        = 360;
            az_vec   = linspace(-pi/2, pi/2, N);
            el_vec   = linspace(-pi/2, pi/2, N);
            [az, el] = meshgrid(az_vec, el_vec);
            m        = obj.TaperExpHorz;
            n        = obj.TaperExpVert;
            Graw     = cos(az/2).^(2*m) .* cos(el/2).^(2*n);
            dAz      = az_vec(2) - az_vec(1);
            dEl      = el_vec(2) - el_vec(1);
            totalPow = sum(Graw .* cos(el), 'all') * dAz * dEl;
            nf       = 4*pi / totalPow;
        end

    end
end
