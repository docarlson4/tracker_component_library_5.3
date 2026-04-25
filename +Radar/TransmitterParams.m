classdef TransmitterParams < handle
% TransmitterParams  Radar transmitter parameters and link-budget quantities.
%
%   Holds all scalar parameters needed by the range equation and waveform
%   builders.  No signal generation.  Analogous to ReceiverNoise.
%
%   CONSTRUCTION
%     p = TransmitterParams(fc, Pt, tau, PRF)
%     p = TransmitterParams(fc, Pt, tau, PRF, Name=Value, ...)
%
%   Required
%     fc   – carrier frequency      [Hz]
%     Pt   – peak transmit power    [W]
%     tau  – pulse width            [s]
%     PRF  – pulse repetition freq. [Hz]
%
%   Optional Name-Value
%     bandwidth  – tx bandwidth override  [Hz]    (default: 1/tau)
%     waveform   – 'rect'|'lfm'|'barker'|'custom'  (default: 'lfm')
%     phaseCode  – chip phase vector [rad]  (required for barker/custom)
%     losses_dB  – ohmic + distribution losses  [dB]  (default: 0)
%
%   EXAMPLE
%     p = TransmitterParams(10e9, 100e3, 10e-6, 1e3, bandwidth=5e6)
%     p.lambda,  p.TBP,  p.Runamb

    properties (SetAccess = private)
        fc          % Carrier frequency       [Hz]
        Pt          % Peak transmit power     [W]
        tau         % Pulse width             [s]
        PRF         % Pulse repetition freq.  [Hz]
        bandwidth   % Instantaneous tx BW     [Hz]
        waveform    % 'rect' | 'lfm' | 'barker' | 'custom'
        phaseCode   % Chip phase vector [rad] (column)
        losses_dB   % Transmitter losses      [dB]
    end

    properties (Dependent)
        lambda                % Wavelength                    [m]
        PRI                   % Pulse repetition interval     [s]
        dutyCycle             % Duty cycle                    [-]
        Pavg                  % Average power                 [W]
        Ep                    % Pulse energy                  [J]
        Runamb                % Unambiguous range             [m]
        fdamb                 % Unambiguous Doppler half-ext. [Hz]
        timeBandwidthProduct  % TBP = bandwidth × tau
        EIRP_dBW              % EIRP (no antenna gain)        [dBW]
    end

    properties (Constant, Access = private)
        c0 = 299792458;
    end

    % ==================================================================
    methods

        function obj = TransmitterParams(fc, Pt, tau, PRF, opts)
            arguments
                fc              (1,1) double {mustBePositive}
                Pt              (1,1) double {mustBePositive}
                tau             (1,1) double {mustBePositive}
                PRF             (1,1) double {mustBePositive}
                opts.bandwidth  (1,1) double {mustBeNonnegative} = 0
                opts.waveform   (1,:) char   {mustBeMember(opts.waveform, ...
                                    {'rect','lfm','barker','custom'})} = 'lfm'
                opts.phaseCode  (:,1) double = pi*[0 0 0 0 0 1 1 0 0 1 0 1 0].'
                opts.losses_dB  (1,1) double {mustBeNonnegative}  = 0
            end
            obj.fc        = fc;
            obj.Pt        = Pt;
            obj.tau       = tau;
            obj.PRF       = PRF;
            obj.waveform  = opts.waveform;
            obj.phaseCode = opts.phaseCode;
            obj.losses_dB = opts.losses_dB;
            obj.bandwidth = opts.bandwidth;
            if obj.bandwidth == 0
                obj.bandwidth = 1 / tau;
            end
            assert(tau < 1/PRF, ...
                'TransmitterParams: tau >= PRI (duty cycle >= 1).');
            assert(~strcmp(obj.waveform,'lfm') || obj.bandwidth*tau > 1+1e-9, ...
                'TransmitterParams: LFM requires TBP > 1.');
            assert(~ismember(obj.waveform,{'barker','custom'}) || ~isempty(obj.phaseCode), ...
                'TransmitterParams: phaseCode required for ''%s''.', obj.waveform);
        end

        % ── Dependent getters ─────────────────────────────────────────
        function v = get.lambda(obj),  v = obj.c0 / obj.fc;           end
        function v = get.PRI(obj),     v = 1 / obj.PRF;               end
        function v = get.dutyCycle(obj),v= obj.tau * obj.PRF;         end
        function v = get.Pavg(obj),    v = obj.Pt * obj.dutyCycle;    end
        function v = get.Ep(obj),      v = obj.Pt * obj.tau;          end
        function v = get.Runamb(obj),  v = obj.c0 / (2*obj.PRF);     end
        function v = get.fdamb(obj),   v = obj.PRF / 2;               end
        function v = get.timeBandwidthProduct(obj), v = obj.bandwidth * obj.tau; end
        function v = get.EIRP_dBW(obj), v = 10*log10(obj.Pt) - obj.losses_dB; end

        function disp(obj)
            fprintf('  TransmitterParams\n');
            fprintf('    fc        = %s\n',       Radar.TransmitterParams.fmtHz(obj.fc));
            fprintf('    lambda    = %.4g mm\n',  obj.lambda*1e3);
            fprintf('    Pt        = %.4g kW\n',  obj.Pt/1e3);
            fprintf('    tau       = %.4g µs\n',  obj.tau*1e6);
            fprintf('    PRF       = %.4g Hz\n',  obj.PRF);
            fprintf('    bandwidth = %s\n',       Radar.TransmitterParams.fmtHz(obj.bandwidth));
            fprintf('    TBP       = %.4g\n',     obj.timeBandwidthProduct);
            fprintf('    duty      = %.3g%%\n',   obj.dutyCycle*100);
            fprintf('    Pavg      = %.4g W\n',   obj.Pavg);
            fprintf('    Ep        = %.4g µJ\n',  obj.Ep*1e6);
            fprintf('    Runamb    = %.4g km\n',  obj.Runamb/1e3);
            fprintf('    ±fdamb    = %.4g Hz\n',  obj.fdamb);
            fprintf('    EIRP      = %.4g dBW\n', obj.EIRP_dBW);
            fprintf('    waveform  = %s\n',       obj.waveform);
            if ismember(obj.waveform,{'barker','custom'})
                fprintf('    code len  = %d chips\n', numel(obj.phaseCode));
            end
            fprintf('    losses    = %.4g dB\n',  obj.losses_dB);
        end

    end

    methods (Static, Access = private)
        function s = fmtHz(f)
            if     f >= 1e9, s = sprintf('%.4g GHz', f/1e9);
            elseif f >= 1e6, s = sprintf('%.4g MHz', f/1e6);
            elseif f >= 1e3, s = sprintf('%.4g kHz', f/1e3);
            else,            s = sprintf('%.4g Hz',  f);
            end
        end
    end

end
