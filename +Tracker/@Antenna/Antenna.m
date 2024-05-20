classdef Antenna
    %ANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ApertureNormal
        AzimuthCosineN
        ElevationCosineN
    end
    
    methods
        function obj = Antenna(varargin)
            %ANTENNA Construct an instance of this class
            % Propagation vector in ANT CS
            %   k_hat = [cos(el).*sin(az); sin(el); cos(el).*cos(az)]

            DEFAULT.ApertureNormal   = [0;0;1]; % ANT CS
            DEFAULT.AzimuthCosineN   = 2;
            DEFAULT.ElevationCosineN = 2;

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end
        end
        
        function G = Gain(obj,az,el)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            m = obj.AzimuthCosineN;
            n = obj.ElevationCosineN;
            f = 4*gamma(m+1)/gamma(m+0.5)*gamma(n+1.5)/gamma(n+1);
            G = f*cos(az).^(2*m).*cos(el).^(2*n);
        end
    end
end

