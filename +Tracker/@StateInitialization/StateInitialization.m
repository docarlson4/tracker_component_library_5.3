classdef StateInitialization < handle
    %STATEINITIALIZATION Handles various forms of state initialization
    
    properties(Access=private)
        expectedTypes = ["One-Point", "Two-Point", "Three-Point"];
        motionModelTypes = ["NCV", "NCA", "GMV"]
        %         stateDim
        kBuf
        zBuf
        SBuf
        tBuf
    end
    properties
        Type        % input ["One-Point", "Two-Point"]
%         SpaceDim    % input
        VelMin
        VelMax
        AccMax
        JrkMax
        MotionModelObject
    end
    
    methods
        function obj = StateInitialization(varargin)
            % STATEINITIALIZATION
            DEFAULT.Type = "One-Point";
            % DEFAULT.SpaceDim = 3;
            DEFAULT.VelMin = 0;
            DEFAULT.VelMax = 100;
            DEFAULT.AccMax = 10;
            DEFAULT.JrkMax = 1;
            DEFAULT.MotionModelObject = "NCV";

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end
            
            obj = init(obj);
            % clear TwoPoint
            % clear ThreePoint
        end

        function [xNew, SNew, lgclNew] = Initialize(obj, tCur, zCur, SRCur)
            switch obj.Type
                case "One-Point"
                    [xNew, SNew, lgclNew] = OnePoint(obj, zCur, SRCur);
                case "Two-Point"
                    [xNew, SNew, lgclNew] = TwoPoint(obj, tCur, zCur, SRCur);
                case "Three-Point"
                    [xNew, SNew, lgclNew] = ThreePoint(obj, tCur, zCur, SRCur);
            end
        end
    end

    methods (Access = private)
        [xNew, SNew, lgclNew] = OnePoint(obj, zCur, SRCur)
        [xNew, SNew, lgclNew] = TwoPoint(obj, tCur, zCur, SRCur)
        [xNew, SNew, lgclNew] = ThreePoint(obj, tCur, zCur, SRCur)
    end
end

function obj = init(obj)
%INIT Validate inputs
if ~any(contains(obj.Type, obj.expectedTypes))
    error("Wrong Type. Expected " + strjoin(obj.expectedTypes, ", "))
end
if ~any(contains(obj.MotionModelObject.Type, obj.motionModelTypes ))
    error("Wrong MotionModelObject. Expected " ...
        + strjoin(obj.motionModelTypes, ", "))
end
% obj.stateDim = obj.MotionModelObject.StateDim;
obj.kBuf = 0;
end
