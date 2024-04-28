classdef MeasurementModel < handle
    %MEASUREMENTMODEL generates a measurement model. The range and
    %range-rate output variables, "Range" and "RangeRate", repectively, are
    %bistatic range and range-rate. Direction cosines, "DirCosU" and
    %DirCosV" are defined in the antenna coordinate system (ACS). The input
    %variables "tgt_state", "xmtr_state", and "rcvr_state" must have
    %minimum dimension (6 X N), where coordinates (1:3 X N) are position
    %and (4:6 X N) are velocity in a global Cartesian coordinate system,
    %e.g. ENU, NED, ECEF, etc. The set of rotation matrices, "M" (3 X 3 X
    %N), transform the direction cosines to the ACS. Care must be taken to
    %distinguish PORT, FORWARD, or STARBOARD pointing antenna systems when
    %defining M.
    %
    % - MeasurementModel uses the following TCL functionality
    %       Cart2Ruv
    %       getRangeRate
    %       rangeGradient
    %       uvGradient
    %       rangeRateGradient
    
    properties
        % INPUT: The 6XN position and velocity vector of the target
        tgt_state
        % INPUT: The 6XN position and velocity vector of the transmitter
        xmtr_state
        % INPUT: The 6XN position and velocity vector of the receiver
        rcvr_state
        % INPUT: 3X3XN rotation matrices from the global C.S. to the receiver
        M
        % INPUT: Any combination of ["Range","DirCosU","DirCosV","RangeRate"]
        types

        z % measurement vector (nz X N)
        R % measurement covariance matrix (nz X nz X N)
        H % state to measurement transformation (nz X ns X N)
        h % measurement function handel h = @(s) H*s
    end
    properties (Access = private)
        nz
        ns
        N
        types_available = ["Range","DirCosU","DirCosV","RangeRate"]
    end

    methods
        % Make sure the type of motion model is avaialble
        function verify_types(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            % Verify input types
            for k = 1:numel(obj.types)
                if ~matches(obj.types_available, obj.types(k))
                    str = sprintf('%s, %s, %s, %s', obj.types_available');
                    error('Wrong type. Input \"%s\" must be one of the \"%s\".', ...
                        obj.types(k), str)
                end
            end

            % Verify Xmtr and Rcvr sizes match the target
            [nsXmtr, NXmtr] = size(obj.xmtr_state);
            [nsRcvr, NRcvr] = size(obj.rcvr_state);
            [nsTgt, NTgt] = size(obj.tgt_state);
            if (NXmtr ~= NRcvr) || (nsXmtr ~= nsRcvr)
                error('Xmtr and Rcvr have differnet sizes')
            end
            if (NXmtr < NTgt)% || 
                error('Too few Xmtr states')
            end
            if (nsXmtr ~= nsTgt)
                error('Xmtr, Rcvr, and Target have differnet sizes')
            end
        end

        function obj = MeasurementModel(sTgt, sXmtr, sRcvr, M, types)
            %MEASUREMENTMODEL Construct an instance of this class
            %   Detailed explanation goes here
            if (nargin < 1 || isempty(sTgt))
                error('A 6XN, N>0, target state vector must be defined')
            end

            [obj.ns, obj.N] = size(sTgt);
            if obj.ns < 6
                error('Target state vector dimension must be at least 6')
            end
            obj.tgt_state = sTgt;

            if (nargin < 5 || isempty(types))
                obj.types = ["Range", "DirCosU", "DirCosV", "RangeRate"];
            else
                obj.types = types;
            end
            obj.nz = numel(obj.types);

            if (nargin < 4 || isempty(M))
            	obj.M = repmat( eye(3,3), [1, 1, obj.N] );
            else
                obj.M = M(:,:,1:obj.N);
            end

            if (nargin < 3 || isempty(sRcvr))
                obj.rcvr_state = zeros(6, obj.N);
            else
                obj.rcvr_state = sRcvr;
            end

            if (nargin < 2 || isempty(sXmtr))
                obj.xmtr_state = zeros(6, obj.N);
            else
                obj.xmtr_state = sXmtr;
            end

            % Make sure the type of motion model is avaialble
            verify_types(obj)

            obj.xmtr_state = obj.xmtr_state(:, 1:obj.N);
            obj.rcvr_state = obj.rcvr_state(:, 1:obj.N);

            useHalfRange = false;
            J=zeros(4, 6, obj.N);

            for kN = 1:obj.N
                J(1, 1:3, kN) = rangeGradient(obj.tgt_state(1:3,kN), ...
                    useHalfRange, ...
                    obj.xmtr_state(1:3,kN), obj.rcvr_state(1:3,kN));

                J(2:3, 1:3, kN) = uvGradient(obj.tgt_state(1:3,kN), ...
                    obj.rcvr_state(1:3,kN),obj.M(:,:,kN));

                J(4, :, kN) = rangeRateGradient(obj.tgt_state(1:6, kN), ...
                    useHalfRange, ...
                    obj.xmtr_state(1:6,kN), obj.rcvr_state(1:6,kN));
            end
            obj.H = J;

            includeW = true;
            zRuv = Cart2Ruv(obj.tgt_state, useHalfRange, ...
                obj.xmtr_state, obj.rcvr_state, obj.M, includeW);
            zRR = getRangeRate(obj.tgt_state(1:6,:), useHalfRange, ...
                obj.xmtr_state, obj.rcvr_state);
            obj.z = [zRuv; zRR];
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

