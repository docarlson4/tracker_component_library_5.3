classdef RadarTransmitter < handle
% RadarTransmitter  Parametric radar transmitter with waveform generation.
%
%   Supports rectangular, LFM (chirp), and phase-coded (Barker/custom)
%   waveforms.  No toolbox dependencies.  Requires R2021a+.
%
%   CONSTRUCTION
%     tx = RadarTransmitter(fc, Pt, tau, PRF)
%     tx = RadarTransmitter(fc, Pt, tau, PRF, Name=Value, ...)
%
%   Required arguments
%     fc   – carrier frequency        [Hz]
%     Pt   – peak transmit power      [W]
%     tau  – pulse width              [s]
%     PRF  – pulse repetition freq.   [Hz]
%
%   Optional Name-Value pairs
%     bandwidth  – override tx bandwidth [Hz]  (default: 1/tau)
%     waveform   – 'rect' | 'lfm' | 'barker' | 'custom'  (default: 'rect')
%     phaseCode  – phase vector [rad] (required for 'barker'/'custom')
%     losses_dB  – ohmic + distribution losses (default: 0)
%
%   EXAMPLE — X-band LFM pulse
%     tx = RadarTransmitter(10e9, 100e3, 10e-6, 1e3, ...
%                           bandwidth=5e6, waveform='lfm');
%     [s, t] = tx.generatePulse(50e6);
%     h      = tx.matchedFilter(50e6);
%     [chi, tau_ax, fd_ax] = tx.ambiguityFunction(50e6);
%
%   EXAMPLE — Barker-13 phase code
%     b13 = pi * [0 0 0 0 0 1 1 0 0 1 0 1 0]';
%     tx  = RadarTransmitter(3e9, 10e3, 26e-6, 500, ...
%                            waveform='barker', phaseCode=b13);

    % ------------------------------------------------------------------
    properties (SetAccess = private)
        fc          % Carrier frequency [Hz]
        Pt          % Peak transmit power [W]
        tau         % Pulse width [s]
        PRF         % Pulse repetition frequency [Hz]
        bandwidth   % Instantaneous tx bandwidth [Hz]
        waveform    % 'rect' | 'lfm' | 'barker' | 'custom'
        phaseCode   % Phase code vector [rad] (column)
        losses_dB   % Transmitter losses [dB]
    end

    % ------------------------------------------------------------------
    properties (Dependent)
        lambda                % Wavelength [m]
        PRI                   % Pulse repetition interval [s]
        dutyCycle             % Duty cycle [dimensionless]
        Pavg                  % Average power [W]
        Ep                    % Pulse energy [J]
        Runamb                % Unambiguous range [m]
        fdamb                 % Unambiguous Doppler half-extent [Hz]  (+/-fdamb)
        timeBandwidthProduct  % TBP = bandwidth x tau
        EIRP_dBW              % EIRP relative to 1 W (no antenna gain)
    end

    % ------------------------------------------------------------------
    properties (Constant, Access = private)
        c0 = 299792458;   % Speed of light [m/s]
    end

    % ==================================================================
    methods

        % --------------------------------------------------------------
        function obj = RadarTransmitter(fc, Pt, tau, PRF, options)
            arguments
                fc          (1,1) double {mustBePositive}
                Pt          (1,1) double {mustBePositive}
                tau         (1,1) double {mustBePositive}
                PRF         (1,1) double {mustBePositive}
                options.bandwidth  (1,1) double {mustBeNonnegative} = 0
                options.waveform   (1,:) char   {mustBeMember(options.waveform, ...
                                            {'rect','lfm','barker','custom'})} = 'rect'
                options.phaseCode  (:,1) double = pi * [0 0 0 0 0 1 1 0 0 1 0 1 0]';
                options.losses_dB  (1,1) double = 0
            end

            obj.fc        = fc;
            obj.Pt        = Pt;
            obj.tau       = tau;
            obj.PRF       = PRF;
            obj.waveform  = options.waveform;
            obj.phaseCode = options.phaseCode;
            obj.losses_dB = options.losses_dB;
            obj.bandwidth = options.bandwidth;
            if obj.bandwidth == 0
                obj.bandwidth = 1 / tau;
            end

            assert(tau < 1/PRF, ...
                'RadarTransmitter: tau >= PRI (duty cycle >= 1).');
            assert(~strcmp(obj.waveform,'lfm') || obj.bandwidth*tau > 1+1e-9, ...
                'RadarTransmitter: LFM requires TBP > 1.');
            assert(~ismember(obj.waveform,{'barker','custom'}) || ~isempty(obj.phaseCode), ...
                'RadarTransmitter: phaseCode required for ''%s''.', obj.waveform);
        end

        % --------------------------------------------------------------
        % Dependent getters
        % --------------------------------------------------------------
        function v = get.lambda(obj),  v = obj.c0 / obj.fc;           end
        function v = get.PRI(obj),     v = 1 / obj.PRF;               end
        function v = get.dutyCycle(obj),v = obj.tau * obj.PRF;        end
        function v = get.Pavg(obj),    v = obj.Pt * obj.dutyCycle;    end
        function v = get.Ep(obj),      v = obj.Pt * obj.tau;          end
        function v = get.Runamb(obj),  v = obj.c0 / (2 * obj.PRF);   end
        function v = get.fdamb(obj),   v = obj.PRF / 2;               end
        function v = get.timeBandwidthProduct(obj), v = obj.bandwidth * obj.tau; end
        function v = get.EIRP_dBW(obj)
            v = 10*log10(obj.Pt) - obj.losses_dB;
        end

        % --------------------------------------------------------------
        function [s, t] = generatePulse(obj, fs)
        % generatePulse(fs)  ->  unit-amplitude complex envelope, sampled at fs.
        %
        %   s : Nx1 complex, |s(t)| = 1 before loss scaling
        %   t : Nx1 time axis [s], t in [0, tau)
        %
        %   Scale s by sqrt(Pt) for volts into a 1-ohm load.
            arguments
                obj
                fs (1,1) double {mustBePositive}
            end

            t  = (0 : 1/fs : obj.tau - 1/fs).';
            N  = numel(t);

            switch obj.waveform
                case 'rect'
                    s = Radar.RadarTransmitter.waveRect(N);
                case 'lfm'
                    s = Radar.RadarTransmitter.waveLFM(t, obj.tau, obj.bandwidth);
                case {'barker','custom'}
                    s = Radar.RadarTransmitter.wavePhaseCode(N, obj.phaseCode);
            end

            s = s * 10^(-obj.losses_dB / 20);
        end

        % --------------------------------------------------------------
        function h = matchedFilter(obj, fs)
        % matchedFilter(fs)  ->  MF impulse response (time-reversed conjugate).
        %
        %   Normalised so the on-axis correlation peak equals 1.
            arguments
                obj
                fs (1,1) double {mustBePositive}
            end

            s = obj.generatePulse(fs);
            h = conj(flipud(s));
            h = h / real(s' * s);
        end

        % --------------------------------------------------------------
        function [chi, delayAxis, dopplerAxis] = ambiguityFunction(obj, fs, nDoppler)
        % ambiguityFunction(fs)  ->  narrowband ambiguity |chi(tau,fd)|, normalised.
        %
        %   chi         : (2N-1) x nDoppler, peak = 1
        %   delayAxis   : (2N-1)x1  [s]
        %   dopplerAxis : 1xnDoppler [Hz], symmetric +/-PRF/2
            arguments
                obj
                fs       (1,1) double {mustBePositive}
                nDoppler (1,1) double {mustBePositive, mustBeInteger} = 256
            end

            nDoppler = nDoppler + mod(nDoppler, 2);

            [s, t] = obj.generatePulse(fs);
            N      = numel(s);

            dopplerAxis = (-nDoppler/2:nDoppler/2-1) * (fs/nDoppler);
            delayAxis   = (-(N-1) : (N-1)).' / fs;

            sSteered = s .* exp(1j * 2*pi * dopplerAxis .* t);

            Nfft = 2^nextpow2(2*N - 1);
            S    = fft(s,        Nfft);
            SS   = fft(sSteered, Nfft, 1);
            corr = ifft(conj(S) .* SS, [], 1);

            chi = [corr(Nfft-N+2:Nfft, :); corr(1:N, :)];
            chi = abs(chi) / max(abs(corr(:)));
        end

        % --------------------------------------------------------------
        function [Pout_dB, f, BW] = spectrum(obj, fs, nFFT)
        % spectrum(fs, nFFT)  ->  PSD [dB, peak=0] and 3-dB BW estimate [Hz].
            arguments
                obj
                fs   (1,1) double {mustBePositive}
                nFFT (1,1) double {mustBeNonnegative, mustBeInteger} = 0
            end

            if nFFT == 0
                nFFT = 2^nextpow2(16 * round(obj.tau * fs));
            end
            s = obj.generatePulse(fs);
            S       = fftshift(fft(s, nFFT));
            Sp      = abs(S).^2;
            Pout_dB = 10*log10(Sp / max(Sp));
            BW      = (sum(Sp >= max(Sp)/2) / nFFT) * fs;
            f = (-nFFT/2:nFFT/2-1)*fs/nFFT;
        end

        % --------------------------------------------------------------
        function disp(obj)
            fprintf('  RadarTransmitter\n');
            fprintf('    fc           = %s\n', obj.fmtFreq(obj.fc));
            fprintf('    lambda       = %.4g mm\n',  obj.lambda * 1e3);
            fprintf('    Pt           = %.4g kW\n',  obj.Pt / 1e3);
            fprintf('    tau          = %.4g µs\n',  obj.tau * 1e6);
            fprintf('    PRF          = %.4g Hz\n',  obj.PRF);
            fprintf('    bandwidth    = %s\n', obj.fmtFreq(obj.bandwidth));
            fprintf('    TBP          = %.4g\n',     obj.timeBandwidthProduct);
            fprintf('    duty cycle   = %.3g%%\n',   obj.dutyCycle * 100);
            fprintf('    Pavg         = %.4g W\n',   obj.Pavg);
            fprintf('    Ep           = %.4g µJ\n',  obj.Ep * 1e6);
            fprintf('    Runamb       = %.4g km\n',  obj.Runamb / 1e3);
            fprintf('    +/-fdamb     = %.4g Hz\n',  obj.fdamb);
            fprintf('    EIRP         = %.4g dBW\n', obj.EIRP_dBW);
            fprintf('    waveform     = %s\n',       obj.waveform);
            if ismember(obj.waveform, {'barker','custom'})
                fprintf('    code length  = %d chips\n', numel(obj.phaseCode));
            end
            fprintf('    losses       = %.4g dB\n',  obj.losses_dB);
        end

    end   % public methods

    % ==================================================================
    methods (Static, Access = private)

        % --------------------------------------------------------------
        function s = waveRect(N)
        % Unit-amplitude rectangular envelope; N samples.
            s = ones(N, 1);
        end

        % --------------------------------------------------------------
        function s = waveLFM(t, tau, B)
        % Linear FM (chirp) envelope.
        %
        %   Instantaneous frequency sweeps linearly over [-B/2, +B/2]
        %   centred on the pulse midpoint.  Quadratic phase argument:
        %     phi(t) = pi * (B/tau) * (t - tau/2)^2
            tc = t - tau/2;
            s  = exp(1j * pi * (B/tau) * tc.^2);
        end

        % --------------------------------------------------------------
        function s = wavePhaseCode(N, code)
        % Phase-coded envelope from a column vector of chip phases [rad].
        %
        %   Chip boundaries are placed at round(k*N/Nc) for k = 0..Nc,
        %   distributing the N samples evenly with at most 1-sample
        %   rounding per boundary — no tail truncation or zero-padding.
            Nc       = numel(code);
            edges    = round((0:Nc) * (N/Nc));   % (Nc+1) boundary indices
            chipLens = diff(edges);               % samples per chip, sum == N
            s        = exp(1j * repelem(code, chipLens));
        end

        % --------------------------------------------------------------
        function str = fmtFreq(f)
            if     f >= 1e9, str = sprintf('%.4g GHz', f/1e9);
            elseif f >= 1e6, str = sprintf('%.4g MHz', f/1e6);
            elseif f >= 1e3, str = sprintf('%.4g kHz', f/1e3);
            else,            str = sprintf('%.4g Hz',  f);
            end
        end

    end

end
