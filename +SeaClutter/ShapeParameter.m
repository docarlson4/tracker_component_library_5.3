classdef ShapeParameter < handle
    properties
        Model (1,1) string = "GIT"
        % Data (:,1) double = []
        GrazingAngle (1,1) double {mustBePositive} = 5 * pi/180;
        ResCellArea (1,1) double {mustBePositive} = 100;
        SwellDir (1,1) double {mustBeNumeric} = 0; % between beam and swell
        Polarization (1,1) string = "VV"
        Enabled (1,1) logical = true
    end

    properties (Dependent, GetAccess=public, SetAccess=private)
        Shape
    end

    properties(Access=private)
        kPol
    end

    methods
        function obj = ShapeParameter(options)
            % Define valid arguments dynamically from class properties
            arguments
                options.Model (1,1) string = "GIT"
                % Data (:,1) double = []
                options.GrazingAngle (1,1) double {mustBePositive} = 5 * pi/180;
                options.ResCellArea (1,1) double {mustBePositive} = 100;
                options.SwellDir (1,1) double {mustBeNumeric} = 0; % between beam and swell
                options.Polarization (1,1) string = "VV"
                options.Enabled (1,1) logical = true
            end

            % Assign all provided name-value pairs
            propNames = fieldnames(options);
            for i = 1:numel(propNames)
                obj.(propNames{i}) = options.(propNames{i});
            end

            switch obj.Polarization
                case 'VV'
                    obj.kPol = 1.39;
                case 'HH'
                    obj.kPol = 2.09;
            end

        end

        function value = get.Shape(obj)
            value = 2*log10(obj.GrazingAngle)/3 ...
                + 5*log10(obj.ResCellArea)/8 ...
                - cos(2*obj.SwellDir)/3 ...
                - obj.kPol;
            value = 10.^(value);
        end
        
    end
end
