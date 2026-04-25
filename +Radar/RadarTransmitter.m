classdef RadarTransmitter < handle
% RadarTransmitter  Facade composing TransmitterParams and TransmitterWaveform.
%
%   Mirrors the RadarReceiver pattern: sub-objects are exposed as public
%   properties for direct access; instance methods are thin delegation wrappers.
%
%   CONSTRUCTION
%     tx = RadarTransmitter(params)
%     tx = RadarTransmitter(fc, Pt, tau, PRF)            – builds params internally
%     tx = RadarTransmitter(fc, Pt, tau, PRF, Name=Value, ...)
%
%   SUB-OBJECT ACCESS
%     tx.params   → TransmitterParams   (fc, Pt, tau, PRF, bandwidth, TBP, …)
%     tx.wave     → TransmitterWaveform (static; access waveform builders directly)
%
%   WAVEFORM METHODS  (convenience wrappers)
%     s          = tx.generatePulse(fs)
%     h          = tx.matchedFilter(fs)
%     [chi,…]    = tx.ambiguityFunction(fs, nDoppler)
%     [Pdb,f,BW] = tx.spectrum(fs, nFFT)
%
%   PASS-THROUGH PROPERTIES  (for range equation compatibility with RadarReceiver)
%     tx.fc, tx.Pt, tx.tau, tx.PRF, tx.bandwidth, tx.lambda,
%     tx.timeBandwidthProduct, tx.losses_dB
%
%   EXAMPLE
%     % --- via facade ---
%     tx = RadarTransmitter(10e9, 100e3, 10e-6, 1e3, bandwidth=5e6, waveform='lfm');
%     disp(tx)
%     [s, t]           = tx.generatePulse(50e6);
%     [chi, tau_ax, fd_ax] = tx.ambiguityFunction(50e6);
%     TransmitterWaveform.plotAmbiguity(chi, tau_ax, fd_ax);
%
%     % --- sub-objects directly ---
%     p  = TransmitterParams(10e9, 100e3, 10e-6, 1e3, bandwidth=5e6);
%     t  = (0 : 1/50e6 : p.tau - 1/50e6).';
%     s  = TransmitterWaveform.lfm(t, p.tau, p.bandwidth);
%     h  = TransmitterWaveform.matchedFilter(s);
%     tx = RadarTransmitter(p);

    % ------------------------------------------------------------------
    properties (SetAccess = private)
        params    % TransmitterParams
        wave      % TransmitterWaveform class handle (for discoverability)
    end

    % ── Pass-through dependent properties (range equation compatibility) ─
    properties (Dependent)
        fc; Pt; tau; PRF; bandwidth; losses_dB
        lambda; timeBandwidthProduct
    end

    % ==================================================================
    methods

        function obj = RadarTransmitter(varargin)
        % RadarTransmitter(params)
        % RadarTransmitter(fc, Pt, tau, PRF, ...)
            if nargin >= 1 && isa(varargin{1}, 'Radar.TransmitterParams')
                obj.params = varargin{1};
            else
                % Forward all arguments to TransmitterParams
                obj.params = Radar.TransmitterParams(varargin{:});
            end
            obj.wave = Radar.TransmitterWaveform;
        end

        % ── Pass-through getters ──────────────────────────────────────
        function v = get.fc(obj),                   v = obj.params.fc;                   end
        function v = get.Pt(obj),                   v = obj.params.Pt;                   end
        function v = get.tau(obj),                  v = obj.params.tau;                  end
        function v = get.PRF(obj),                  v = obj.params.PRF;                  end
        function v = get.bandwidth(obj),            v = obj.params.bandwidth;            end
        function v = get.losses_dB(obj),            v = obj.params.losses_dB;            end
        function v = get.lambda(obj),               v = obj.params.lambda;               end
        function v = get.timeBandwidthProduct(obj), v = obj.params.timeBandwidthProduct; end

        % ── Waveform convenience wrappers ─────────────────────────────

        function [s, t] = generatePulse(obj, fs)
        % generatePulse(fs)  ->  [s, t]  complex envelope + time axis at rate fs.
            arguments
                obj
                fs (1,1) double {mustBePositive}
            end
            t = (0 : 1/fs : obj.params.tau - 1/fs).';
            s = Radar.TransmitterWaveform.pulse(obj.params, fs);
        end

        function h = matchedFilter(obj, fs)
        % matchedFilter(fs)  ->  MF impulse response.
            arguments
                obj
                fs (1,1) double {mustBePositive}
            end
            s = Radar.TransmitterWaveform.pulse(obj.params, fs);
            h = Radar.TransmitterWaveform.matchedFilter(s);
        end

        function [chi, delayAxis, dopplerAxis] = ambiguityFunction(obj, fs, nDoppler)
        % ambiguityFunction(fs, nDoppler)  ->  [chi, delayAxis, dopplerAxis].
            arguments
                obj
                fs       (1,1) double {mustBePositive}
                nDoppler (1,1) double {mustBePositive, mustBeInteger} = 256
            end
            s = Radar.TransmitterWaveform.pulse(obj.params, fs);
            [chi, delayAxis, dopplerAxis] = ...
                Radar.TransmitterWaveform.ambiguityFunction(s, fs, nDoppler);
        end

        function [Pdb, f, BW] = spectrum(obj, fs, nFFT)
        % spectrum(fs, nFFT)  ->  [Pdb, f, BW].
            arguments
                obj
                fs   (1,1) double {mustBePositive}
                nFFT (1,1) double {mustBeNonnegative, mustBeInteger} = 0
            end
            s = Radar.TransmitterWaveform.pulse(obj.params, fs);
            [Pdb, f, BW] = Radar.TransmitterWaveform.spectrum(s, fs, nFFT);
        end

        % ── Display ───────────────────────────────────────────────────

        function disp(obj)
            fprintf('  RadarTransmitter\n');
            disp(obj.params);
        end

    end

end
