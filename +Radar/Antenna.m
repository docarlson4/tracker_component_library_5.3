classdef Antenna < handle
    %ANTENNA  Phased-array antenna model.
    %
    %   Separates the element pattern from the array geometry through
    %   composition: any AntennaElement subclass can be plugged in.
    %
    %   CONSTRUCTION
    %     a = Antenna()
    %     a = Antenna(Name=Value, ...)
    %
    %   NAMED PARAMETERS (with defaults)
    %     Element       AntennaElement   – element pattern model          [CosineElement()]
    %     Wavelength    (1,1) double > 0 – operating wavelength, m        [0.03]
    %     DeltaElemHorz (1,1) double > 0 – horizontal element spacing, m  [λ/2]
    %     DeltaElemVert (1,1) double > 0 – vertical element spacing, m    [λ/2]
    %     NumSubHorz    (1,1) int    > 0 – subarrays horizontally         [4]
    %     NumSubVert    (1,1) int    > 0 – subarrays vertically           [1]
    %     NumElemHorz   (1,1) int    > 0 – elements per subarray (horiz.) [2]
    %     NumElemVert   (1,1) int    > 0 – elements per subarray (vert.)  [8]
    %     AzScan        (1,1) double     – azimuth scan angle, rad        [0]
    %     ElScan        (1,1) double     – elevation scan angle, rad      [0]
    %
    %   ELEMENT TYPES
    %     CosineElement        – default; cos^(2N) separable model
    %     PatchElement         – rectangular microstrip patch
    %     VivaldiElement       – tapered-slot (Vivaldi) approximation
    %     CrossedDipoleElement – two orthogonal dipoles in quadrature
    %     NumericalElement     – data-driven pattern via interpolation
    %
    %   EXAMPLE
    %     % X-band Vivaldi array, 8×8 elements, λ/2 spacing
    %     a = Antenna( ...
    %           Element      = VivaldiElement(TaperExpHorz=1.5), ...
    %           Wavelength   = 0.03, ...
    %           NumElemHorz  = 8, ...
    %           NumElemVert  = 8, ...
    %           NumSubHorz   = 1, ...
    %           NumSubVert   = 1, ...
    %           AzScan       = 30*pi/180 );
    %     a.PlotFullPat();
    %
    %   Propagation vector convention (antenna coordinate system):
    %     k_hat = [cos(el)*sin(az); sin(el); cos(el)*cos(az)]

    % -----------------------------------------------------------------
    %  Public properties
    % -----------------------------------------------------------------
    properties
        Element       {mustBeA(Element, 'Radar.Elements.AntennaElement')} = Radar.Elements.CosineElement()
        Wavelength    (1,1) double {mustBePositive}                      = 0.03
        DeltaElemHorz (1,1) double {mustBePositive}                      = 0.015  % λ/2 default
        DeltaElemVert (1,1) double {mustBePositive}                      = 0.015  % λ/2 default
        NumSubHorz    (1,1) double {mustBePositive, mustBeInteger}       = 4
        NumSubVert    (1,1) double {mustBePositive, mustBeInteger}       = 1
        NumElemHorz   (1,1) double {mustBePositive, mustBeInteger}       = 2
        NumElemVert   (1,1) double {mustBePositive, mustBeInteger}       = 8
        AzScan        (1,1) double                                       = 0
        ElScan        (1,1) double                                       = 0
    end

    % -----------------------------------------------------------------
    %  Dependent property — recomputed from Wavelength automatically
    % -----------------------------------------------------------------
    properties (Dependent, SetAccess = private)
        k   % wavenumber 2π/λ (rad/m)
    end

    % =================================================================
    methods
        % =================================================================

        function obj = Antenna(options)
            %ANTENNA  Construct a phased-array Antenna object.
            %
            %   a = Antenna()  – all defaults (CosineElement, λ=0.03 m, etc.)
            %   a = Antenna(Element=PatchElement(), Wavelength=0.05)
            %
            %   Note: DeltaElemHorz/Vert default to 0.015 m = λ/2 at the
            %   default wavelength.  Pass them explicitly when changing λ, or
            %   update them after construction.
            arguments
                options.Element  {mustBeA(options.Element, 'Radar.Elements.AntennaElement')} ...
                    = Radar.Elements.CosineElement()
                options.Wavelength    (1,1) double {mustBePositive}           = 0.03
                options.DeltaElemHorz (1,1) double {mustBePositive}           = 0.015
                options.DeltaElemVert (1,1) double {mustBePositive}           = 0.015
                options.NumSubHorz    (1,1) double {mustBePositive, mustBeInteger} = 4
                options.NumSubVert    (1,1) double {mustBePositive, mustBeInteger} = 1
                options.NumElemHorz   (1,1) double {mustBePositive, mustBeInteger} = 2
                options.NumElemVert   (1,1) double {mustBePositive, mustBeInteger} = 8
                options.AzScan        (1,1) double                            = 0
                options.ElScan        (1,1) double                            = 0
            end
            fields = fieldnames(options);
            for i = 1:numel(fields)
                obj.(fields{i}) = options.(fields{i});
            end
        end

        % ---- Dependent getter --------------------------------------- %
        function k = get.k(obj)
            k = 2*pi / obj.Wavelength;
        end

        % ==============================================================
        %  Gain methods
        % ==============================================================

        function G = ElementGain(obj, az, el)
            %ELEMENTGAIN  Delegate to the attached AntennaElement object.
            %
            %   This thin wrapper preserves backward compatibility: callers
            %   can address ElementGain on the array without knowing which
            %   element type is installed.
            G = obj.Element.ElementGain(az, el);
        end

        function G = SubarrayGain(obj, az, el)
            %SUBARRAYGAIN  ElementGain × subarray array factor.
            dx          = obj.DeltaElemHorz;
            dy          = obj.DeltaElemVert;
            [kx, ky]    = obj.spatialFreqs(az, el);
            af          = obj.arrayFactor(kx, ky, obj.NumElemHorz, obj.NumElemVert, dx, dy);
            G           = obj.ElementGain(az, el) .* abs(af).^2 ...
                / (obj.NumElemHorz * obj.NumElemVert);
        end

        function G = PhasedArrayGain(obj, az, el)
            %PHASEDARRAY GAIN  SubarrayGain × phased-array array factor.
            %
            %   The cos(az)·cos(el) factor is a projected-aperture correction
            %   for the element orientation relative to the look direction.
            dx          = obj.DeltaElemHorz * obj.NumElemHorz;  % subarray pitch
            dy          = obj.DeltaElemVert * obj.NumElemVert;
            [kx, ky]    = obj.spatialFreqs(az, el);
            af          = obj.arrayFactor(kx, ky, obj.NumSubHorz, obj.NumSubVert, dx, dy);
            G           = obj.SubarrayGain(az, el) .* abs(af).^2 ...
                .* cos(az) .* cos(el) ...
                / (obj.NumSubHorz * obj.NumSubVert);
        end

        % ==============================================================
        %  Beamwidth
        % ==============================================================

        function bw = AzBeamwidth(obj)
            %AZBEAMWIDTH  3-dB azimuth beamwidth (radians).
            N  = obj.NumElemHorz * obj.NumSubHorz;
            x  = obj.calcBeamwidth(N);
            bw = 2 * asin(x / (obj.k * obj.DeltaElemHorz));
        end

        function bw = ElBeamwidth(obj)
            %ELBEAMWIDTH  3-dB elevation beamwidth (radians).
            N  = obj.NumElemVert * obj.NumSubVert;
            x  = obj.calcBeamwidth(N);
            bw = 2 * asin(x / (obj.k * obj.DeltaElemVert));
        end

        % ==============================================================
        %  Plots
        % ==============================================================

        function ax = PlotElemPat(obj, az0, el0)
            %PLOTELEMPATH  Element pattern plotted over an az/el grid (degrees).
            arguments
                obj
                az0 (1,:) double = -89:89
                el0 (1,:) double = -89:89
            end
            [~, ~, gain] = obj.evalPattern(@obj.ElementGain, az0, el0, -50);
            ax = obj.renderSurf(az0, el0, gain, ...
                '\bf Azimuth (deg)', '\bf Elevation (deg)', ...
                sprintf('Element Pattern – %s', obj.Element.Description));
        end

        function ax = PlotSubarrayPat(obj, az0, el0)
            %PLOTSUBARRAYPAT  Subarray pattern on u-v (direction-cosine) axes.
            arguments
                obj
                az0 (1,:) double = -89:89
                el0 (1,:) double = -89:89
            end
            [az, el, gain] = obj.evalPattern(@obj.SubarrayGain, az0, el0, -50);
            u  = cos(el) .* sin(az);
            v  = sin(el);
            ax = obj.renderSurf(u, v, gain, '\bf u', '\bf v', 'Subarray Pattern');
        end

        function ax = PlotFullPat(obj, az0, el0)
            %PLOTFULLPAT  Full array pattern on u-v (direction-cosine) axes.
            arguments
                obj
                az0 (1,:) double = -89:0.1:89
                el0 (1,:) double = -89:0.1:89
            end
            [az, el, gain] = obj.evalPattern(@obj.PhasedArrayGain, az0, el0, 0);
            u  = cos(el) .* sin(az);
            v  = sin(el);
            ax = obj.renderSurf(u, v, gain, '\bf u', '\bf v', 'Full Array Pattern');
        end

    end % public methods

    % =================================================================
    methods (Access = private)
        % =================================================================

        function [kx, ky] = spatialFreqs(obj, az, el)
            %SPATIALFREQS  Steering-corrected spatial frequencies.
            kx = obj.k .* (sin(az).*cos(el)  - sin(obj.AzScan).*cos(obj.ElScan));
            ky = obj.k .* (sin(el)            - sin(obj.ElScan));
        end

        function af = arrayFactor(~, kx, ky, Nx, Ny, dx, dy)
            %ARRAYFACTOR  2-D uniform rectangular array factor, centred at origin.
            idxX = (1:Nx) - mean(1:Nx);
            idxY = (1:Ny) - mean(1:Ny);
            af   = zeros(size(kx));
            for cx = idxX
                for cy = idxY
                    af = af + exp(1i*(kx*cx*dx + ky*cy*dy));
                end
            end
        end

        function x = calcBeamwidth(~, N)
            %CALCBEAMWIDTH  Newton-Raphson: sinc_array(x,N) = 1/√2.
            f  = @(x) sin(N*x/2) ./ (N*sin(x/2)) - 1/sqrt(2);
            fp = @(x) (N*sin(x/2).*cos(N*x/2) - cos(x/2).*sin(N*x/2)) ...
                ./ (2*N*sin(x/2).^2);
            x  = 0.1;
            for iter = 1:100
                dx = f(x) ./ fp(x);
                x  = x - dx;
                if abs(dx) < 1e-6, break; end
            end
        end

        function [az, el, gain] = evalPattern(obj, gainFcn, az0, el0, floorDB)
            %EVALPATTERN  Meshgrid + gain evaluation + dB conversion.
            deg      = pi / 180;
            [az, el] = meshgrid(az0*deg, el0*deg);
            gain     = dB(gainFcn(az, el), floorDB);
        end

        function ax = renderSurf(~, x, y, gain, xLabel, yLabel, titleStr)
            %RENDERSURF  Shared surface-plot formatter.
            figure;
            surf(x, y, gain);
            shading flat;
            xlabel(xLabel);
            ylabel(yLabel);
            maxG = max(gain(:));
            title({ sprintf('\\bf\\fontsize{14} %s (dB)', titleStr), ...
                sprintf('\\bf\\fontsize{12} Max Gain = %.1f dB', maxG) });
            ax = gca;
        end

    end % private methods

end % classdef
