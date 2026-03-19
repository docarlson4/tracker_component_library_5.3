classdef CrossedDipoleElement < Radar.Elements.AntennaElement
%CROSSEDDIPOLEELEMENT  Two orthogonal half-wave dipoles fed in phase quadrature.
%
%   PATTERN MODEL
%
%   Dipole 1 lies along the array horizontal axis (x).  Its power pattern
%   is proportional to 1 − (cos(el)·sin(az))², the squared projection of
%   the observation unit vector onto x.
%
%   Dipole 2 lies along the array vertical axis (y).  Its pattern is
%   proportional to 1 − sin²(el) = cos²(el).
%
%   For equal-amplitude, 90°-quadrature feeding (RHCP or LHCP), the total
%   radiated power pattern is the incoherent sum of both dipoles:
%
%     G(az,el) = 0.75 · [1 + cos²(el)·cos²(az)]
%
%   (The factor 0.75 = 1.5/2 comes from the isotropic-normalised half-wave
%   dipole gain of 1.64 → 1.5 on average, split equally across two ports.)
%   Peak gain at broadside is 1.5 (1.76 dBi), consistent with a single
%   half-wave dipole — the second dipole fills the orthogonal null.
%
%   CONSTRUCTION
%     e = CrossedDipoleElement()
%
%   Reference: Balanis, "Antenna Theory," 3rd ed., §4.6 (crossed dipoles).

    methods

        function obj = CrossedDipoleElement()
        %CROSSEDDIPOLEELEMENT  No configurable parameters.
            arguments; end
            obj.Description = "Crossed dipole (RHCP/LHCP quadrature feed)";
        end

        function G = ElementGain(~, az, el)
        %ELEMENTGAIN  G = 0.75*(1 + cos²(el)*cos²(az))
        %
        %   Peak = 1.5 (1.76 dBi) at broadside.
        %   Pattern is hemispherically symmetric — no ground-plane model.
            G = 0.75 .* (1 + cos(el).^2 .* cos(az).^2);
        end

    end
end
