classdef TargetSwarm < handle
    %TARGETSWARM Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        expectedTypes = [
            "Type1"
            "Type2"
            ]
        km = 1e3
        mph = 0.44704
        deg = pi/180
        ddm2rrm = diag([pi/180, pi/180, 1])
        pos_plat_geo
        pos_plat_ecf
        pos_ref_geo
        pos_ref_ecf
        EC
        CE
    end

    properties
        Type
        RadarObject
        RadarLocationGEO
        ScanNumber
        RangeLimits
        VelocityLimits
        TargetNumber
        SwarmReferenceGEO
    end

    methods
        function obj = TargetSwarm(varargin)
            %TARGETSWARM Construct an instance of this class
            obj.Type = "Type1";
            obj.RadarObject = Tracker.RadarReceiver;
            obj.RadarLocationGEO = [35;-85;300];
            obj.SwarmReferenceGEO = [35;-85;300];
            obj.ScanNumber = 50;
            obj.RangeLimits = [1, 50] * obj.km;
            obj.VelocityLimits = [200, 600] * obj.mph;
            obj.TargetNumber = 3;

            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, obj.(fn) ), ...
                fieldnames( obj ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end

            % Make sure the type of motion model is avaialble
            verify_type(obj)

            % Initialization
            obj = init(obj);

            % Generate swarm
            obj = genSwarm(obj);

        end

        function ShowTypes(obj)
            disp(obj.expectedTypes)
        end
    end

    methods (Access = private)
        % Generate swarm
        obj = genSwarm(obj);

        % Initialization
        function obj = init(obj)
            % Reference location
            obj.pos_ref_geo = obj.ddm2rrm * obj.SwarmReferenceGEO;
            obj.pos_ref_ecf = ellips2Cart(obj.pos_ref_geo);
            % Radar platform location
            obj.pos_plat_geo = obj.ddm2rrm * obj.RadarLocationGEO;
            obj.pos_plat_ecf = ellips2Cart(obj.pos_plat_geo);
            obj.EC = getENUAxes(obj.pos_plat_geo);
            obj.CE = obj.EC';
        end
        
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
    end
end
% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-11-03 12:31
