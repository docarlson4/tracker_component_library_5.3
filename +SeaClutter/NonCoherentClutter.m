classdef NonCoherentClutter < handle
    %NONCOHERENTCLUTTER This first level corresponds to the generation of
    %sea clutter data intended for non-Doppler processing using a single
    %channel antenna: thus, the CPI is equal to one single PRI (NFFT = 1).
    %
    % CONSTRUCTION
    %   a = NonCoherentClutter()
    %   a = NonCoherentClutter(Name=Value, ...)
    % 
    % NAMED PARAMETERS (with defaults)
    %   TextureModel (1,1) string   - ["Gamma", "InverseGamma"]
    %   ShapeParam (1,1) double     - [nu]
    %   CorrCoeffRng (1,1) function_handle = @deal
    % 
    % 
    % 
    % 
    % 
    % [0] Kemkemian, S., Lupinski, L., Corretja, V., Cottron, R., & Watts,
    % S. (2015, May). Performance assessment of multi-channel radars using
    % simulated sea clutter. In 2015 IEEE Radar Conference (RadarCon) (pp.
    % 1015-1020). IEEE.

    % -----------------------------------------------------------------
    %  Public properties
    % -----------------------------------------------------------------
    properties
        TextureModel %(1,1) string = "Gamma" % ["Gamma", "InverseGamma"]
        ShapeParam %(1,1) double = 1    % [nu]
        CorrCoeffRng %(1,1) function_handle = @deal
        NumRange %(1,1) double {mustBePositive} = 64;
        NumAzimuth %(1,1) double {mustBePositive} = 128;
    end

    % -----------------------------------------------------------------
    %  Private properties
    % -----------------------------------------------------------------
    properties (Access = private)
    end

    % -----------------------------------------------------------------
    %  Dependent property — recomputed from Wavelength automatically
    % -----------------------------------------------------------------
    properties (Dependent, GetAccess = public, SetAccess = private)
        ScaleParam
    end

    methods
        function obj = NonCoherentClutter(opts)
            %NONCOHERENTCLUTTER Construct an instance of this class
            arguments
                opts.TextureModel (1,1) string = "Gamma" % ["Gamma", "InverseGamma"]
                opts.ShapeParam (1,1) double {mustBePositive} = 1.0 % [nu]
                opts.CorrCoeffRng (1,1) function_handle = @deal
                opts.NumRange (1,1) double {mustBePositive} = 64;
                opts.NumAzimuth (1,1) double {mustBePositive} = 128;
            end
            fields = fieldnames(opts);
            for i = 1:numel(fields)
                obj.(fields{i}) = opts.(fields{i});
            end

        end

        % ---- Dependent getter --------------------------------------- %
        function value = get.ScaleParam(obj)
            value = 1/obj.ShapeParam;
        end

    end % public methods

    methods (Access = private)
    end % private methods

end % classdef

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2026-03-20 13:40
