classdef PatchElement < Radar.Elements.AntennaElement
%PATCHELEMENT  Rectangular microstrip patch antenna element pattern.
%
%   Uses the two-slot transmission-line model.  The patch is oriented with
%   its resonant dimension (L) along the array horizontal axis and its
%   non-resonant dimension (W) along the array vertical axis.
%
%   PATTERN MODEL (approximate, ignores substrate effects)
%
%     E-plane  (az-scan, el=0): F_E(az) = cos(π·Lr·sin(az))
%     H-plane  (el-scan, az=0): F_H(el) = sinc(π·Wr·sin(el))·cos(el)
%
%   where Lr = L/λ, Wr = W/λ, and MATLAB's sinc(x) = sin(πx)/(πx).
%
%   The 2-D pattern combines both factors with a cos(el) half-space taper:
%     G(az,el) ∝ cos²(π·Lr·sin(az)·cos(el)) · sinc²(π·Wr·sin(el)) · cos²(el)
%
%   Normalisation is computed numerically at construction.
%
%   CONSTRUCTION
%     e = PatchElement()
%     e = PatchElement(PatchLengthRatio=0.47, PatchWidthRatio=0.6)
%
%   PROPERTIES
%     PatchLengthRatio  L/λ, resonant dimension (default 0.49 ≈ λ/2)
%     PatchWidthRatio   W/λ, non-resonant dimension (default 0.6)

    properties
        PatchLengthRatio (1,1) double {mustBePositive} = 0.49
        PatchWidthRatio  (1,1) double {mustBePositive} = 0.60
    end

    properties (Access = private)
        NormFactor (1,1) double = 1   % numerical normalisation constant
    end

    methods

        function obj = PatchElement(options)
            arguments
                options.PatchLengthRatio (1,1) double {mustBePositive} = 0.49
                options.PatchWidthRatio  (1,1) double {mustBePositive} = 0.60
            end
            obj.PatchLengthRatio = options.PatchLengthRatio;
            obj.PatchWidthRatio  = options.PatchWidthRatio;
            obj.NormFactor       = obj.computeNorm();
            obj.Description      = sprintf( ...
                "Patch element (L/λ=%.4g, W/λ=%.4g)", ...
                obj.PatchLengthRatio, obj.PatchWidthRatio);
        end

        function G = ElementGain(obj, az, el)
        %ELEMENTGAIN  Patch pattern normalised to integrate to 4π sr.
        %
        %   Returns zero for el outside (-π/2, π/2) to model the ground
        %   plane suppressing back-lobe radiation.
            Lr = obj.PatchLengthRatio;
            Wr = obj.PatchWidthRatio;

            % E-plane factor: cos of electrical-length projection
            F_E = cos(pi .* Lr .* sin(az) .* cos(el));

            % H-plane factor: sinc of slot-width projection, with cos(el) taper
            F_H = sinc(Wr .* sin(el)) .* cos(el);

            G = obj.NormFactor .* F_E.^2 .* F_H.^2;

            % Suppress back hemisphere (ground plane)
            G(el < -pi/2 | el > pi/2) = 0;
            G(az < -pi/2 | az > pi/2) = 0;
        end

    end

    % -----------------------------------------------------------------
    methods (Access = private)
    % -----------------------------------------------------------------

        function nf = computeNorm(obj)
        %COMPUTENORM  Integrate unnormalised pattern over the front hemisphere
        %   and return the scale factor that makes total radiated power = 4π.
            N        = 360;
            az_vec   = linspace(-pi/2, pi/2, N);
            el_vec   = linspace(-pi/2, pi/2, N);
            [az, el] = meshgrid(az_vec, el_vec);
            Lr       = obj.PatchLengthRatio;
            Wr       = obj.PatchWidthRatio;
            F_E      = cos(pi .* Lr .* sin(az) .* cos(el));
            F_H      = sinc(Wr .* sin(el)) .* cos(el);
            Graw     = F_E.^2 .* F_H.^2;
            % Spherical-coordinate area element: dΩ = cos(el) daz del
            dAz = az_vec(2) - az_vec(1);
            dEl = el_vec(2) - el_vec(1);
            totalPow = sum(Graw .* cos(el), 'all') * dAz * dEl;
            nf = 4*pi / totalPow;
        end

    end
end
