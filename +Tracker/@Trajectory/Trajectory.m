classdef Trajectory < handle & matlab.mixin.Heterogeneous
    %TRAJECTORY defines the state vector "s" of a moving object in an
    %East-North-Up (ENU) coordinate system.

    properties
        T           % INPUT: revisit time
        v0          % INPUT: constant speed along legs and turns
        pos0        % INPUT: inital position
        leg_dirs    % INPUT: leg directions of compass rose
        leg_times   % INPUT: time spent on each leg
        turn_accs   % INPUT: turn acceleration, in g's, between legs

        Ns          % OUTPUT: number of discrete states
        s           % OUTPUT: s(:,1) = [pos; vel; ang speed], (7 X Ns)
        mode        % OUTPUT: {+1,0,-1} == {PORT, LINEAR, STARBOARD} motion
        t           % OUTPUT: time domain (1 X Ns) in seconds
        leg_angs    % OUTPUT: body frame yaw angles
    end

    properties (Access = private)
        %% Compass directions and angles in True North CS
        compass_dir = ["N","NNE","NE","ENE","E","ESE","SE","SSE", ...
            "S","SSW","SW","WSW","W","WNW","NW","NNW"];
        compass_ang = (0:22.5:359);
    end

    methods
        function obj = Trajectory(varargin)
            %TRAJECTORY Construct an instance of this class
            %   Detailed explanation goes here
            p = inputParser;

            addParameter(p,'revisit_time',1,@isnumeric);
            addParameter(p,'speed',20,@isnumeric);
            addParameter(p,'initial_position',[0;0;0],@isnumeric);
            check_legs = @(x) all(ismember(x,obj.compass_dir));
            addParameter(p,'leg_directions',["N", "NE"],check_legs);
            addParameter(p,'leg_times',[60,45],@isnumeric);
            addParameter(p,'turn_accelerations',0.1,@isnumeric);

            parse(p,varargin{:})

            obj.T = p.Results.revisit_time;
            obj.v0 = p.Results.speed;
            obj.pos0 = p.Results.initial_position;
            obj.leg_dirs = p.Results.leg_directions;
            obj.leg_times = p.Results.leg_times;
            obj.turn_accs = p.Results.turn_accelerations * 9.81;

            validate_entries(obj);

            obj = gen_traj(obj);
        end
    end

    methods(Access = private)
        function obj = gen_traj(obj)
            deg = pi/180;

            num_leg = numel(obj.leg_dirs);
            num_seg = 2*num_leg-1;

            modes = false(1,num_seg); modes(2:2:num_seg) = true;
            times = nan(1,num_seg); times(1:2:num_seg) = obj.leg_times;
            %acc = zeros(1,num_seg); acc(2:2:num_seg) = obj.turn_accs;

            ang = zeros(1,num_leg);
            for kL = 1:num_leg
                ang(kL) = obj.compass_ang(matches(obj.compass_dir, obj.leg_dirs(kL)));
            end
            del_phi = -diff(ang)*deg;

            for kA = 1:num_leg-1

                if obj.turn_accs(kA) > 0 && del_phi(kA) < 0
                    del_phi(kA) = del_phi(kA) + 2*pi;
                elseif obj.turn_accs(kA) < 0 && del_phi(kA) > 0
                    del_phi(kA) = 2*pi - del_phi(kA);
                end
            end
            
            N_turn = ceil( abs(del_phi)*obj.v0./(abs(obj.turn_accs)*obj.T) );
            N = times/obj.T; N(2:2:num_seg) = N_turn;
            %times(2:2:num_seg) = N_turn*T;
            w = zeros(1,num_seg); w(2:2:num_seg) = sign(obj.turn_accs).*abs(del_phi)./(obj.T*N_turn);
            N = ceil(N);
            num_scan = sum(N)+1;
            N_cum = [0, cumsum(N)]+1;
            pos = obj.pos0;
            heading = [sin(ang(1)*deg);cos(ang(1)*deg);0];
            vel = obj.v0*heading;

            obj.s = [pos;vel;0]*ones(1,num_scan);
            obj.t = zeros(1,num_scan);
            obj.mode = zeros(1,num_scan);
            idx = [1:2, 4:5, 7];
            for kS = 2:num_seg+1

                for k = N_cum(kS-1)+1:N_cum(kS)
                    obj.s(7,k-1) = w(kS-1);
                    obj.s(idx,k) = FCoordTurn2D(obj.T,obj.s(idx,k-1)) * obj.s(idx,k-1);
                    obj.t(k) = (k-1)*obj.T;
                    %mode(k) = modes(kS-1);
                    obj.mode(k) = sign(w(kS-1));
                end
            end

            % Gather leg_angs
            lgcl_angs = false(1,numel(obj.compass_ang));
            for k = 1:numel(obj.leg_dirs)
                lgcl_angs = lgcl_angs | ismember(obj.compass_dir, obj.leg_dirs(k));
            end
            obj.leg_angs = obj.compass_ang(lgcl_angs);
            
            obj.Ns = numel(obj.t);
%             if center
%                 obj.s(1:2,:) = obj.s(1:2,:) - (max(obj.s(1:2,:),[],2) + min(obj.s(1:2,:),[],2))/2;
%             end
        end

        function validate_entries(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            num_legs = numel(obj.leg_dirs);
            if num_legs ~= numel(obj.leg_times)
                error('The number of \"leg_times\" must equal the number of \"leg_directions\"')
            end
            if num_legs-1 ~= numel(obj.turn_accs)
                error('The number of \"turn_accelerations\" must be one less than the the number of \"leg_directions\"')
            end
        end
    end
end

