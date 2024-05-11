classdef StateInitialization < handle
    %STATEINITIALIZATION Handles various forms of state initialization
    
    properties(Access=private)
        expectedTypes = ["One-Point", "Two-Point", "Three-Point"];
        motionModelTypes = ["NCV", "NCA", "GMV"]
    end
    properties
        Type        % input ["One-Point", "Two-Point"]
        StateDim    % input xDim
        SpaceDim    % input
        VelMin
        VelMax
        AccMax
        JrkMax
        MotionModelType
    end
    
    methods
        function obj = StateInitialization(varargin)
            % STATEINITIALIZATION
            DEFAULT.Type = "One-Point";
            DEFAULT.StateDim = 6;
            DEFAULT.SpaceDim = 3;
            DEFAULT.VelMin = 0;
            DEFAULT.VelMax = 100;
            DEFAULT.AccMax = 10;
            DEFAULT.JrkMax = 1;
            DEFAULT.MotionModelType = "NCV";

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end
            
            obj = init(obj);
            clear TwoPoint
            clear ThreePoint
        end

        function [xNew, SNew] = Initialize(obj, tCur, zCur, SRCur)
            switch obj.Type
                case "One-Point"
                    [xNew, SNew] = OnePoint(obj, zCur, SRCur);
                case "Two-Point"
                    [xNew, SNew] = TwoPoint(obj, tCur, zCur, SRCur);
                case "Three-Point"
                    [xNew, SNew] = ThreePoint(obj, tCur, zCur, SRCur);
            end
        end
    end

    methods (Access = private)
        [xNew, SNew] = OnePoint(obj, zCur, SRCur)
        [xNew, SNew] = TwoPoint(obj, tCur, zCur, SRCur)
        [xNew, SNew] = ThreePoint(obj, tCur, zCur, SRCur)
    end
end

function obj = init(obj)
%INIT Validate inputs
if ~any(contains(obj.Type, obj.expectedTypes))
    error("Wrong Type. Expected " + strjoin(obj.expectedTypes, ", "))
end
if ~any(contains(obj.MotionModelType, obj.motionModelTypes ))
    error("Wrong MotionModelType. Expected " ...
        + strjoin(obj.motionModelTypes, ", "))
end
end
