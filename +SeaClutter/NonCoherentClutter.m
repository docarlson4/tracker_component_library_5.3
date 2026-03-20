classdef NonCoherentClutter < handle
    %NONCOHERENTCLUTTER This first level corresponds to the generation of
    %sea clutter data intended for non-Doppler processing using a single
    %channel antenna: thus, the CPI is equal to one single PRI (NFFT = 1).
    %
    %   CONSTRUCTION
    %
    %   NAMED PARAMETERS (with defaults)
    %
    % [0] Kemkemian, S., Lupinski, L., Corretja, V., Cottron, R., & Watts,
    % S. (2015, May). Performance assessment of multi-channel radars using
    % simulated sea clutter. In 2015 IEEE Radar Conference (RadarCon) (pp.
    % 1015-1020). IEEE.

    % -----------------------------------------------------------------
    %  Public properties
    % -----------------------------------------------------------------
    properties
        Type
        Property1
    end

    % -----------------------------------------------------------------
    %  Private properties
    % -----------------------------------------------------------------
    properties (Access = private)
        expectedTypes = [
            "Type1"
            "Type2"
            ]
    end

    % -----------------------------------------------------------------
    %  Dependent property — recomputed from Wavelength automatically
    % -----------------------------------------------------------------
    properties (Dependent, SetAccess = private)
        k   % wavenumber 2π/λ (rad/m)
    end

    methods
        function obj = NonCoherentClutter(options)
            %NONCOHERENTCLUTTER Construct an instance of this class
            arguments
                options.Type = "Type1";
                options.Property1 = "Property1";
            end
            fields = fieldnames(options);
            for i = 1:numel(fields)
                obj.(fields{i}) = options.(fields{i});
            end

            % Make sure the type of motion model is avaialble
            verify_type(obj)
        end

        % ---- Dependent getter --------------------------------------- %
        function k = get.k(obj)
            k = 2*pi / obj.Wavelength;
        end

        function ShowTypes(obj)
            disp(obj.expectedTypes)
        end
    end % public methods

    methods (Access = private)
        % Make sure the type of motion model is available
        function verify_type(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if ~contains(obj.expectedTypes, obj.Type)
                strJ = strjoin(obj.expectedTypes,", ");
                str = sprintf('%s', strJ);
                error('Wrong type. Input Type \"%s\" must be one of \"%s\".', ...
                    obj.Type, str)
            end
        end
    end % private methods

end % classdef

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2026-03-20 13:40
