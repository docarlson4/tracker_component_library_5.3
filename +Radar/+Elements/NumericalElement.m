classdef NumericalElement < Radar.Elements.AntennaElement
%NUMERICALELEMENT  Antenna element pattern defined by measured or simulated data.
%
%   Stores a 2-D gain table on a regular az/el grid and uses bilinear
%   interpolation (griddedInterpolant) for queries at arbitrary angles.
%   Gain is clamped to zero outside the supplied grid extents.
%
%   CONSTRUCTION
%     e = NumericalElement(azGrid, elGrid, gainData)
%
%   INPUT
%     azGrid   (1,Na) double – azimuth sample points, radians, monotonic
%     elGrid   (1,Ne) double – elevation sample points, radians, monotonic
%     gainData (Ne,Na) double – linear power gain; rows = elevation,
%                                cols = azimuth.  Must be non-negative.
%
%   OPTIONS (name-value)
%     Description  string   – human-readable label (default: 'Numerical')
%     Interpolation string  – griddedInterpolant method (default: 'linear')
%
%   STATIC FACTORY
%     e = NumericalElement.fromFile(filename)
%       Loads a .mat file containing variables 'azGrid', 'elGrid', and
%       'gainData' (same conventions as above).
%
%   EXAMPLE
%     az  = linspace(-pi/2, pi/2, 181);
%     el  = linspace(-pi/2, pi/2, 91);
%     [A,E] = meshgrid(az, el);
%     G     = 4 .* cos(A).^2 .* cos(E).^2;   % cosine-squared stand-in
%     elem  = NumericalElement(az, el, G, Description="My sim data");

    properties (SetAccess = private)
        AzGrid   (1,:) double   % azimuth sample grid (rad)
        ElGrid   (1,:) double   % elevation sample grid (rad)
        GainData (:,:) double   % linear gain table [Ne × Na]
    end

    properties (Access = private)
        Interpolant   % griddedInterpolant object (built at construction)
    end

    % -----------------------------------------------------------------
    methods
    % -----------------------------------------------------------------

        function obj = NumericalElement(azGrid, elGrid, gainData, options)
        %NUMERICALELEMENT  Construct from az/el grid and gain table.
            arguments
                azGrid   (1,:) double {mustBeVector, mustBeFinite}
                elGrid   (1,:) double {mustBeVector, mustBeFinite}
                gainData (:,:) double {mustBeNonnegative}
                options.Description  (1,1) string  = "Numerical element"
                options.Interpolation (1,1) string = "linear"
            end

            % Validate grid/data size consistency
            assert(size(gainData,1) == numel(elGrid) && ...
                   size(gainData,2) == numel(azGrid), ...
                'NumericalElement: gainData must be [numel(elGrid) × numel(azGrid)].');

            % Ensure grids are strictly monotonic (required by griddedInterpolant)
            assert(all(diff(azGrid) > 0), 'azGrid must be strictly increasing.');
            assert(all(diff(elGrid) > 0), 'elGrid must be strictly increasing.');

            obj.AzGrid      = azGrid;
            obj.ElGrid      = elGrid;
            obj.GainData    = gainData;
            obj.Description = options.Description;
            obj.Interpolant = griddedInterpolant( ...
                {elGrid, azGrid}, gainData, ...
                options.Interpolation, 'none');   % 'none' → NaN outside grid
        end

        function G = ElementGain(obj, az, el)
        %ELEMENTGAIN  Interpolated gain at arbitrary (az, el) in radians.
        %
        %   Query points outside the stored grid return 0 (no extrapolation).
            G = obj.Interpolant(el, az);
            G(isnan(G)) = 0;   % clamp out-of-range to zero
            G = max(G, 0);     % guard against interpolation undershoot
        end

    end

    % -----------------------------------------------------------------
    methods (Static)
    % -----------------------------------------------------------------

        function obj = fromFile(filename, options)
        %FROMFILE  Load a NumericalElement from a .mat file.
        %
        %   The file must contain the variables:
        %     azGrid   – (1×Na) azimuth grid, radians
        %     elGrid   – (1×Ne) elevation grid, radians
        %     gainData – (Ne×Na) linear power gain
        %
        %   Optionally: 'description' string.
        %
        %   USAGE
        %     e = NumericalElement.fromFile('my_pattern.mat')
        %     e = NumericalElement.fromFile('my_pattern.mat', ...
        %                                    Interpolation='makima')
            arguments
                filename (1,1) string {mustBeFile}
                options.Interpolation (1,1) string = "linear"
            end

            S = load(filename, 'azGrid', 'elGrid', 'gainData');

            requiredVars = {'azGrid','elGrid','gainData'};
            missing = requiredVars(~isfield(S, requiredVars));
            assert(isempty(missing), ...
                'NumericalElement.fromFile: file is missing variable(s): %s', ...
                strjoin(missing, ', '));

            desc = filename;   % use filename as default description
            if isfield(S, 'description'), desc = string(S.description); end

            obj = NumericalElement(S.azGrid, S.elGrid, S.gainData, ...
                Description=desc, Interpolation=options.Interpolation);
        end

    end
end
