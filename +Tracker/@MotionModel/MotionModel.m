classdef MotionModel < handle
    %MOTIONMODEL generates a motion model
    % - Uses the following TCL functionality
    %       FPolyKal
    %       QPolyKalDirectDisc
    % - By default, the process noise is zero. To reset it, use the
    %   SetProcessNoise(arg) member function. The argument is the maximum
    %   change of the highest order moment of the state vector. E.g. for
    %   the NCV model, it is the max change in velocity.
    %
    properties (Access = private)
        expectedTypes = ["NCV","NCA","2DCT","GMV"]
        % sigmaV2 = 0;
    end
    properties
        Type        % input ["NCV","NCA","2DCT","GMV"]
        StateDim    % input xDim
        SpaceDim    % input
        RevisitTime % input sample time
        MaxKinVar   % input NCV, GMV: maxAcc, NCA: maxJrk
        Tau         % input GMV: decorrelation time (s)

        Q  % state process noise covariance (ns X ns)
        SQ % Cholesky decomposition of Q (ns X ns)
        F  % state propagator matrix (ns X ns)
        f  % state propagator function handel f = @(s) F*s
    end

    methods
        function obj = MotionModel(varargin)
            %MOTIONMODEL Construct an instance of this class
            DEFAULT.Type = "NCV";
            DEFAULT.StateDim = 6;
            DEFAULT.SpaceDim = 3;
            DEFAULT.RevisitTime = 1;
            DEFAULT.MaxKinVar = 9.8; % 1 g
            DEFAULT.Tau = 10*DEFAULT.RevisitTime; % s
            
            p = inputParser;

            cellfun( @(fn) p.addOptional( fn, DEFAULT.(fn) ), ...
                fieldnames( DEFAULT ) );

            p.parse(varargin{:});

            param_names = fieldnames( p.Results );
            for k = 1:numel(param_names)
                obj.(param_names{k}) = p.Results.(param_names{k});
            end
            
            % Make sure the type of motion model is avaialble
            verify_type(obj)

            % Set the state vector, process noise, and propagator according
            % to type
            obj = update_motion_model(obj);
        end

        function obj = SetProcessNoise(obj, max_change_highest_order_moment)
            %SETPROCESSNOISE Reset the process noise covariance matrix
            % NCV: mxhom ~ m/s (max change in velocity)
            % NCA: mxhom ~ m/s^2 (max change in acceleration)
            mxhom = max_change_highest_order_moment;
            obj.sigmaV2 = (mxhom/obj.RevisitTime)^2;
            obj = update_motion_model(obj);
        end

    end

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
    
        % Set the state vector, process noise, and propagator according
        % to type
        function obj = update_motion_model(obj)
            switch obj.Type
                case "NCV"
                    order = 1;
                    % Process noise
                    q = processNoiseSuggest('PolyKal-ROT',obj.MaxKinVar, ...
                        obj.RevisitTime);
                    obj.Q = QPolyKalDirectDisc( ...
                        obj.RevisitTime, ...
                        obj.StateDim, ...
                        order, ...
                        q);
                    % Propagator
                    obj.F = FPolyKal(obj.RevisitTime,obj.StateDim,order);
                    obj.f = @(x) obj.F*x;
                case "NCA"
                    order = 2;
                    % Process noise
                    q = processNoiseSuggest('PolyKal-ROT',obj.MaxKinVar, ...
                        obj.RevisitTime);
                    obj.Q = QPolyKalDirectDisc( ...
                        obj.RevisitTime, ...
                        obj.StateDim, ...
                        order, ...
                        q);
                    % Propagator
                    obj.F = FPolyKal(obj.RevisitTime,obj.StateDim,order);
                    obj.f = @(x) obj.F*x;
                case "GMV"
                    %     case "Gauss-Markov Velocity (GMV)"
                    %Parameters for the dynamic model. We are using a first-order Gauss-Markov
                    %model.
                    order = 1;
                    %Rule-of-thumb process noise suggestion.
                    q = processNoiseSuggest('PolyKal-ROT',obj.MaxKinVar, ...
                        obj.RevisitTime);
                    obj.Q = QGaussMarkov( ...
                        obj.RevisitTime, ...
                        obj.StateDim, ...
                        q, ...
                        obj.Tau, ...
                        order);%Process noise covariance matrix.

                    obj.F = FGaussMarkov( ...
                        obj.RevisitTime, ...
                        obj.StateDim, ...
                        obj.Tau, ...
                        order);%State transition matrix

            end
            obj.SQ = cholSemiDef(obj.Q, 'lower'); 
        end
    end
end

