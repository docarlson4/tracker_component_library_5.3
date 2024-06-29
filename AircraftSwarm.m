classdef AircraftSwarm < handle
    %AIRCRAFTSWARM Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        expectedTypes = [
            "FixedWing"
            "Rotary"
            ]
    end

    properties
        Type
        NumAircraft
        NumSamples
        MeanLifetime    % how long target is detected
        MeanFirstDet    % when target is detected
        FrameTime
        RangeMax        % target undetected for rng > RangeMax
        RangeMin        % target undetected for rng < RangeMin
        SpeedMax
        SpeedMin

        ZCart % [pos;vel;w]
        ZPol  % [rng;az;rngRate]
    end

    properties(Dependent)
        
    end

    methods
        function obj = AircraftSwarm(varargin)
            %AIRCRAFTSWARM Construct an instance of this class
            DEFAULT.Type = "FixedWing";
            DEFAULT.NumAircraft = 5;
            DEFAULT.NumSamples = 100;
            DEFAULT.MeanLifetime = 75; % samples
            DEFAULT.MeanFirstDet = 10; % samples
            DEFAULT.FrameTime = 2; % sec
            DEFAULT.RangeMax = 40e3; % m
            DEFAULT.RangeMin = 1e3; % m
            DEFAULT.SpeedMax = 20; % m/s (44.7 mph)
            DEFAULT.SpeedMin = 270; % m/s (604 mph)

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

            % Initalize dependent quantities
            obj = init(obj);

        end
    end

    methods (Access = private)
        % Initalize dependent quantities
        function obj = init(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            Ts = obj.FrameTime;
            f = @(x) fTransCoordTurn2D(Ts,x);

            pois = PoissonD;

            num_tgt = obj.NumAircraft;
            Ns = obj.NumSamples;

            obj.ZCart = cell(Ns,1);
            obj.ZPol = cell(Ns,1);

            % Polar location of targets in NED
            % rng_tgt = Rmin + (Rmax-Rmin)*rand(num_tgt, 1);
            rng_tgt = obj.RangeMax*ones(num_tgt, 1);
            alp_tgt = 2*pi*rand(num_tgt, 1);
            pos = [rng_tgt.*sin(alp_tgt), rng_tgt.*cos(alp_tgt)]';

            % Velocity of targets
            speed = obj.SpeedMin + (obj.SpeedMax-obj.SpeedMin)*rand(num_tgt,1);
            bet_tgt = alp_tgt + (-pi/4 + pi/2*rand(num_tgt,1)) + pi;
            vel = [speed.*sin(bet_tgt), speed.*cos(bet_tgt)]';

            % Angular speed
            rv = cross([pos;zeros(1,num_tgt)],[vel;zeros(1,num_tgt)],1);
            w = rv(3,:)./(rng_tgt.*rng_tgt)';
            % w = w .* (2 * (2*rand(1,num_tgt)-0));
            % rCtr = speed./w';
            % max_acc = rCtr.*w'.^2;

            % Target lifetime
            num_lt = 3 + pois.rand([num_tgt,1],obj.MeanLifetime);
            xCart = cell(num_tgt, 1);
            xPol = cell(num_tgt, 1);
            xDot = cell(num_tgt, 1);
            x0 = [pos;vel;w];
            for kT = 1:num_tgt
                xCart{kT} = x0(:,kT);
                xPol{kT} = [
                    vecnorm( xCart{kT}(1:2,1) )
                    atan2( xCart{kT}(1,1), xCart{kT}(2,1) )
                    getRangeRate(xCart{kT}(1:4,1),true)];
                for kLT = 2:min( num_lt(kT), obj.NumSamples )
                    xCart{kT}(:,kLT) = f( xCart{kT}(:,kLT-1) );
                    xDot{kT}(:,kLT) = (xCart{kT}(:,kLT) - xCart{kT}(:,kLT-1))/Ts;

                    xPol{kT} = [xPol{kT}, [
                        vecnorm( xCart{kT}(1:2,kLT) )
                        atan2( xCart{kT}(1,kLT), xCart{kT}(2,kLT) )
                        getRangeRate(xCart{kT}(1:4,kLT),true)] ];
                end
            end

            % First detection
            num_fd = zeros(num_tgt, 1);
            for kFD = 1:num_tgt
                num_fd(kFD) = 1 + pois.rand(1,obj.MeanFirstDet);
            end

            for kS = 1:Ns
                for kT = 1:num_tgt
                    if num_fd(kT) == kS
                        obj.ZCart{kS} = xCart{kT}(:,1);
                        obj.ZPol{kS} = xPol{kT}(:,1);
                    end
                    if (num_fd(kT) < kS) && (kS <= num_fd(kT) + num_lt(kT))
                        idx = kS - num_fd(kT);
                        if xPol{kT}(1,idx) > obj.RangeMax ...
                                || xPol{kT}(1,idx) < obj.RangeMin
                            continue
                        end
                        obj.ZCart{kS} = [obj.ZCart{kS}, xCart{kT}(:,idx)];
                        obj.ZPol{kS} = [obj.ZPol{kS}, xPol{kT}(:,idx)];
                    end
                end
            end
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-06-29 15:14
