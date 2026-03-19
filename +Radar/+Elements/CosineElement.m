classdef CosineElement < Radar.Elements.AntennaElement
%COSINEELEMENT  Separable cosine^(2N) element-pattern model.
%
%   Models the element gain as:
%     G(az,el) = f(m,n) * cos(az)^(2m) * cos(el)^(2n)
%
%   where f(m,n) is a normalisation constant chosen so that integrating G
%   over the unit sphere yields 4π (i.e., unity radiation efficiency).
%   This is the standard analytical model for aperture-type elements.
%
%   CONSTRUCTION
%     e = CosineElement()
%     e = CosineElement(CosineNHorz=3, CosineNVert=4)

    properties
        CosineNHorz (1,1) double {mustBePositive} = 2  % horizontal exponent (m)
        CosineNVert (1,1) double {mustBePositive} = 2  % vertical exponent   (n)
    end

    methods

        function obj = CosineElement(options)
            arguments
                options.CosineNHorz (1,1) double {mustBePositive} = 2
                options.CosineNVert (1,1) double {mustBePositive} = 2
            end
            obj.CosineNHorz  = options.CosineNHorz;
            obj.CosineNVert  = options.CosineNVert;
            obj.Description  = sprintf("Cosine-N element (m=%.4g, n=%.4g)", ...
                                   obj.CosineNHorz, obj.CosineNVert);
        end

        function G = ElementGain(obj, az, el)
        %ELEMENTGAIN  cos(az)^(2m) * cos(el)^(2n), normalised to 4π sr.
            m = obj.CosineNHorz;
            n = obj.CosineNVert;
            % Normalisation via gamma-function identity for cos^(2k) integrals
            f = 4 * gamma(m+1)/gamma(m+0.5) * gamma(n+1.5)/gamma(n+1);
            G = f .* cos(az).^(2*m) .* cos(el).^(2*n);
        end

    end
end
