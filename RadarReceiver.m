classdef RadarReceiver < handle
    %RADARRECEIVER Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = private)
        expectedTypes = [
            "BlairATF"  % Use of Range-Rate Measurements in ATF
            "Type2"
            ]
    end

    properties(Access=private)
        c0 = 299792458;
        deg = pi/180;
        MHz = 1e6;
        GHz = 1e9;
        kHz = 1e3;
        us = 1e-6;
    end

    properties
        Type
        CenterFreq
        Bandwidth
        Beamwidth
        NumPulse
        FrameTime
    end

    properties(Dependent)
        % Co-dependent
        SNR
        SNRdB
        PRF
        PRI
        % Dependent
        RangeAmb
        RangeRateAmb
        RangeUnc
        AzimuthUnc
        RangeRateUnc
        ProbCorrectUnfolding
    end

    methods
        function obj = RadarReceiver(varargin)
            %RADARRECEIVER Construct an instance of this class
            DEFAULT.Type = "BlairATF";
            DEFAULT.CenterFreq = 3*obj.GHz;
            DEFAULT.Bandwidth = 10*obj.MHz;
            DEFAULT.Beamwidth = 2*obj.deg;
            DEFAULT.NumPulse = 20;
            DEFAULT.SNRdB = 10;
            DEFAULT.PRI = 300*obj.us;
            DEFAULT.FrameTime = 2;

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
    end

    %% Plot Methods
    methods
        function PlotUnfoldingProb(obj)
            disp_name = "T_k = " + num2str(obj.FrameTime) + " s";
            figure
            hold on, grid on, box on
            plot(obj.SNRdB, obj.ProbCorrectUnfolding,"s-", ...
                DisplayName=disp_name)
            ylim([0,1])
            xlabel("\bf SNR (dB) ")
            ylabel("\bf Probability of Correct Unfolding ")
            legend(Location="southeast")
        end
    end
    
    %% Get/Set
    properties (Access = private)
        PRIValue % Single property to store either PRI or PRF
        SNRValue % Property to store SNR in linear form
    end
       
    methods
        %% Unambiguous Range
        function value = get.RangeAmb(obj)
            value = obj.c0/(2*obj.PRF);
        end

        %% Unambiguous Range Rate
        function value = get.RangeRateAmb(obj)
            value = obj.c0/(2*obj.CenterFreq*obj.PRI);
        end

        %% RangeUnc
        function value = get.RangeUnc(obj)
            value = obj.c0./(2*obj.Bandwidth*sqrt(obj.SNR));
        end

        %% AzimuthUnc
        function value = get.AzimuthUnc(obj)
            value = sqrt(pi/4) * obj.Beamwidth./sqrt(obj.SNR);
        end
        
        %% RangeRateUnc
        function value = get.RangeRateUnc(obj)
            value = obj.c0./(2*obj.CenterFreq*obj.NumPulse*obj.PRI*sqrt(obj.SNR));
        end
        
        %% ProbCorrectUnfolding
        function value = get.ProbCorrectUnfolding(obj)
            gam = obj.RangeRateUnc./obj.RangeUnc;
            arg = obj.SNR/4 * obj.NumPulse^2./(1 + (2./(gam*obj.FrameTime)).^2);
            value = ChiSquareD.CDF(arg,1);
        end

        %% PRI
        function value = get.PRI(obj)
            value = obj.PRIValue;
        end
        
        function set.PRI(obj, value)
            if value <= 0
                error('PRI must be positive.');
            end
            obj.PRIValue = value;
        end

        %% PRF
        function value = get.PRF(obj)
            value = 1 ./ obj.PRIValue;
        end
        
        function set.PRF(obj, value)
            if value <= 0
                error('PRF must be positive.');
            end
            obj.PRIValue = 1 ./ value;
        end

        %% SNRLinear
        function value = get.SNR(obj)
            value = obj.SNRValue;
        end
        
        function set.SNR(obj, value)
            if value < 0
                error('SNR must be non-negative.');
            end
            obj.SNRValue = value;
        end

        %% SNRdB
        function value = get.SNRdB(obj)
            value = 10 * log10(obj.SNRValue);
        end
        
        function set.SNRdB(obj, value)
            obj.SNRValue = 10.^(value / 10);
        end
    end
    
end
% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-06-28 22:35
