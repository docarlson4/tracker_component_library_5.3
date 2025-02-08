classdef Antenna < handle
    %ANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Cosine N element pattern
        CosineNVert % vertiacl exponent
        CosineNHorz % horizontal exponent
        % Array pattern
        Wavelength % operational wavelength
        DeltaElem % element spacing
        NumSubVert % number of subarrays vertically
        NumSubHorz % number of subarrays horizontally
        NumElemVert % number of elements per subarray vertically
        NumElemHorz % number of elements per subarray horizontally
    end
    
    methods
        function obj = Antenna(varargin)
            %ANTENNA Construct an instance of this class
            % Propagation vector in ANT CS
            %   k_hat = [cos(el).*sin(az); sin(el); cos(el).*cos(az)]

            DEFAULT.CosineNHorz = 2;
            DEFAULT.CosineNVert = 2;
            DEFAULT.Wavelength = 0.03;
            DEFAULT.DeltaElem = DEFAULT.Wavelength/2;
            DEFAULT.NumSubVert = 1;
            DEFAULT.NumSubHorz = 4;
            DEFAULT.NumElemVert = 8;
            DEFAULT.NumElemHorz = 2;

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end
        end
        
        function G = ElementGain(obj,az,el)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            m = obj.CosineNHorz;
            n = obj.CosineNVert;
            f = 4*gamma(m+1)/gamma(m+0.5)*gamma(n+1.5)/gamma(n+1);
            G = f*cos(az).^(2*m).*cos(el).^(2*n);
        end

        function G = SubarrayGain(obj,az,el)
            dx = obj.DeltaElem;
            dy = obj.DeltaElem;
            k = 2*pi/obj.Wavelength;
            kx = k*sin(az).*cos(el);
            ky = k*sin(el);
            af = zeros(size(az));
            idxC = (1:obj.NumElemHorz);
            idxC = idxC - mean(idxC);
            idxR = (1:obj.NumElemVert);
            idxR = idxR - mean(idxR);
            for kC = idxC
                for kR = idxR
                    af = af + exp(1i*(kx*kC*dx + ky*kR*dy));
                end
            end
            G = obj.ElementGain(az,el) .* af.*conj(af) ...
                / (obj.NumElemHorz*obj.NumElemVert);
        end

        function G = PhasedArrayGain(obj,az,el)
            dx = obj.DeltaElem * obj.NumElemHorz;
            dy = obj.DeltaElem * obj.NumElemVert;
            k = 2*pi/obj.Wavelength;
            kx = k*sin(az).*cos(el);
            ky = k*sin(el);
            af = zeros(size(az));
            idxC = (1:obj.NumSubHorz);
            idxC = idxC - mean(idxC);
            idxR = (1:obj.NumSubVert);
            idxR = idxR - mean(idxR);
            for kC = idxC
                for kR = idxR
                    af = af + exp(1i*(kx*kC*dx + ky*kR*dy));
                end
            end
            G = obj.SubarrayGain(az,el) .* af.*conj(af) ...
            .* cos(az).*cos(el) ...
            / (obj.NumSubHorz*obj.NumSubVert);
        end

        function PlotElemPat(obj,az0,el0)
            if nargin == 1
                az0 = -89:89;
                el0 = -89:89;
            end
            deg = pi/180;
            [az,el] = meshgrid(az0*deg,el0*deg);
            gain = dB( obj.ElementGain(az,el), -50 );
            figure
            surf(az0,el0,gain)
            shading flat
            xlabel('\bf Azimuth (deg) ')
            ylabel('\bf Elevation (deg) ')
            title('\bf\fontsize{14} Element Pattern (dB) ')
        end

        function PlotSubarrayPat(obj,az0,el0)
            if nargin == 1
                az0 = -89:89;
                el0 = -89:89;
            end
            deg = pi/180;
            [az,el] = meshgrid(az0*deg,el0*deg);
            gain = dB( obj.SubarrayGain(az,el), -50 );
            figure
            surf(az0,el0,gain)
            surf(cos(el).*sin(az),sin(el),gain)
            shading flat
            xlabel('\bf Azimuth (deg) ')
            ylabel('\bf Elevation (deg) ')
            title(['\bf\fontsize{14} Subarray Pattern (dB) '; ...
                "\bf\fontsize{12} Max Gain "+num2str(max(max(gain)))+" (dB) "])
        end

        function PlotFullPat(obj,az0,el0)
            if nargin == 1
                az0 = -89:0.1:89;
                el0 = -89:0.1:89;
            end
            deg = pi/180;
            [az,el] = meshgrid(az0*deg,el0*deg);
            gain = dB( obj.PhasedArrayGain(az,el), 0 );
            figure
            surf(az0,el0,gain)
            surf(cos(el).*sin(az),sin(el),gain)
            shading flat
            xlabel('\bf Azimuth (deg) ')
            ylabel('\bf Elevation (deg) ')
            title(['\bf\fontsize{14} Full Array Pattern (dB) '; ...
                "\bf\fontsize{12} Max Gain "+num2str(max(max(gain)))+" (dB) "])
        end
    end
end

