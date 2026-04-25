classdef TransmitterWaveform
% TransmitterWaveform  Waveform generation and analysis utilities.
%
%   Analogous to ReceiverDetection: all methods are static so waveforms can
%   be built directly in scripts without touching a transmitter object.
%   The facade RadarTransmitter wraps these calls for convenience.
%
%   STATIC BUILDERS  (call directly for external signal construction)
%     s = TransmitterWaveform.rect(N)
%     s = TransmitterWaveform.lfm(t, tau, B)
%     s = TransmitterWaveform.phaseCode(N, code)
%     s = TransmitterWaveform.pulse(p, fs)       – dispatches on p.waveform
%
%   STATIC ANALYSIS
%     h                    = TransmitterWaveform.matchedFilter(s)
%     [chi, tau_ax, fd_ax] = TransmitterWaveform.ambiguityFunction(s, fs, nDoppler)
%     [Pdb, f, BW]         = TransmitterWaveform.spectrum(s, fs, nFFT)
%
%   EXAMPLE — build and compress an LFM pulse without a transmitter object
%     B = 10e6;  tau = 20e-6;  fs = 50e6;
%     t = (0 : 1/fs : tau - 1/fs).';
%     s = TransmitterWaveform.lfm(t, tau, B);
%     h = TransmitterWaveform.matchedFilter(s);
%     y = conv(s, h, 'same');
%
%   EXAMPLE — ambiguity function
%     [chi, tau_ax, fd_ax] = TransmitterWaveform.ambiguityFunction(s, fs, 256);

    methods (Static)

        % ── Waveform builders ─────────────────────────────────────────

        function s = rect(N)
        % rect(N)  ->  Nx1 unit-amplitude rectangular envelope.
            s = ones(N, 1);
        end

        function s = lfm(t, tau, B)
        % lfm(t, tau, B)  ->  complex LFM envelope.
        %
        %   Instantaneous frequency sweeps over [-B/2, +B/2] symmetric
        %   about the pulse midpoint.
        %     phi(t) = pi * (B/tau) * (t - tau/2)^2
        %
        %   t   : Nx1 time vector [s]
        %   tau : pulse width [s]
        %   B   : chirp bandwidth [Hz]
            arguments
                t   (:,1) double
                tau (1,1) double {mustBePositive}
                B   (1,1) double {mustBePositive}
            end
            s = exp(1j * pi * (B/tau) * (t - tau/2).^2);
        end

        function s = phaseCode(N, code)
        % phaseCode(N, code)  ->  Nx1 phase-coded envelope.
        %
        %   N    : total sample count
        %   code : Nc×1 chip phases [rad]
        %
        %   Chip boundaries at round(k*N/Nc), distributing N samples evenly.
            arguments
                N    (1,1) double {mustBePositive, mustBeInteger}
                code (:,1) double
            end
            Nc       = numel(code);
            edges    = round((0:Nc) * (N/Nc));
            chipLens = diff(edges);
            s        = exp(1j * repelem(code, chipLens));
        end

        function s = pulse(p, fs)
        % pulse(p, fs)  ->  complex envelope dispatched on p.waveform.
        %
        %   p  : TransmitterParams object
        %   fs : sample rate [Hz]
        %
        %   Returns unit-amplitude (before loss scaling) Nx1 complex envelope.
        %   Applies p.losses_dB as a voltage scaling.
            arguments
                p  Radar.TransmitterParams
                fs (1,1) double {mustBePositive}
            end
            t = (0 : 1/fs : p.tau - 1/fs).';
            N = numel(t);
            switch p.waveform
                case 'rect'
                    s = Radar.TransmitterWaveform.rect(N);
                case 'lfm'
                    s = Radar.TransmitterWaveform.lfm(t, p.tau, p.bandwidth);
                case {'barker','custom'}
                    s = Radar.TransmitterWaveform.phaseCode(N, p.phaseCode);
            end
            s = s * 10^(-p.losses_dB / 20);
        end

        % ── Analysis ──────────────────────────────────────────────────

        function h = matchedFilter(s)
        % matchedFilter(s)  ->  MF impulse response (time-reversed conjugate).
        %
        %   Normalised so the on-axis correlation peak equals 1.
        %   s : Nx1 complex envelope (output of any builder above).
            arguments
                s (:,1) {mustBeNumeric}
            end
            h = conj(flipud(s));
            h = h / real(s' * s);
        end

        function [chi, delayAxis, dopplerAxis] = ambiguityFunction(s, fs, nDoppler)
        % ambiguityFunction(s, fs, nDoppler)
        %
        %   Narrowband ambiguity function |chi(tau,fd)| normalised to unit peak.
        %
        %   s        : Nx1 complex envelope
        %   fs       : sample rate [Hz]
        %   nDoppler : number of Doppler bins (default 256, forced even)
        %
        %   Returns
        %     chi         : (2N-1) × nDoppler, peak = 1
        %     delayAxis   : (2N-1)×1  [s]
        %     dopplerAxis : 1×nDoppler [Hz]
            arguments
                s        (:,1) {mustBeNumeric}
                fs       (1,1) double {mustBePositive}
                nDoppler (1,1) double {mustBePositive, mustBeInteger} = 256
            end
            nDoppler = nDoppler + mod(nDoppler,2);
            N        = numel(s);
            t        = (0:N-1).' / fs;

            dopplerAxis = (-nDoppler/2 : nDoppler/2-1) * (fs/nDoppler);
            delayAxis   = (-(N-1) : (N-1)).' / fs;

            sSteered = s .* exp(1j*2*pi * dopplerAxis .* t);

            Nfft = 2^nextpow2(2*N - 1);
            S    = fft(s,        Nfft);
            SS   = fft(sSteered, Nfft, 1);
            corr = ifft(conj(S) .* SS, [], 1);

            chi = [corr(Nfft-N+2:Nfft,:); corr(1:N,:)];
            chi = abs(chi) / max(abs(corr(:)));
        end

        function [Pdb, f, BW] = spectrum(s, fs, nFFT)
        % spectrum(s, fs, nFFT)  ->  two-sided PSD [dB, peak=0] and 3-dB BW [Hz].
        %
        %   nFFT defaults to 2^nextpow2(16*N).
            arguments
                s    (:,1) {mustBeNumeric}
                fs   (1,1) double {mustBePositive}
                nFFT (1,1) double {mustBeNonnegative, mustBeInteger} = 0
            end
            if nFFT == 0
                nFFT = 2^nextpow2(16*numel(s));
            end
            S    = fftshift(fft(s, nFFT));
            Sp   = abs(S).^2;
            Pdb  = 10*log10(Sp / max(Sp));
            BW   = (sum(Sp >= max(Sp)/2) / nFFT) * fs;
            f    = (-nFFT/2 : nFFT/2-1) * (fs/nFFT);
        end

        % ── Plotting ──────────────────────────────────────────────────

        function plotAmbiguity(chi, delayAxis, dopplerAxis, titleStr)
        % plotAmbiguity(chi, delayAxis, dopplerAxis)  — surface + delay/Doppler cuts.
            arguments
                chi         (:,:) double
                delayAxis   (:,1) double
                dopplerAxis (1,:) double
                titleStr    (1,:) char = ''
            end
            chi_dB = 20*log10(chi + eps);

            figure('Name','Ambiguity Function','Position',[60 60 1200 500]);
            tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
            if ~isempty(titleStr)
                title(tl, titleStr, 'FontSize', 11);
            end

            col = [0.2 0.5 0.9];

            nexttile;
            imagesc(dopplerAxis/1e3, delayAxis*1e6, chi_dB);
            clim([-60 0]); colorbar; colormap('jet');
            xlabel('Doppler [kHz]'); ylabel('Delay [µs]');
            title('|χ(τ,f_d)| [dB]');

            nexttile;
            [~, iz] = min(abs(dopplerAxis));
            plot(delayAxis*1e6, chi_dB(:,iz), 'Color',col,'LineWidth',1.2);
            xlabel('Delay [µs]'); ylabel('dB'); grid on;
            title('Range cut  (f_d = 0)');
            ylim([-60 5]);

            nexttile;
            [~, iz] = min(abs(delayAxis));
            plot(dopplerAxis/1e3, chi_dB(iz,:), 'Color',col,'LineWidth',1.2);
            xlabel('Doppler [kHz]'); ylabel('dB'); grid on;
            title('Doppler cut  (τ = 0)');
            ylim([-60 5]);
        end

        function plotSpectrum(Pdb, f, B, titleStr)
        % plotSpectrum(Pdb, f, B)  — annotated PSD plot.
            arguments
                Pdb     (:,1) double
                f       (:,1) double
                B       (1,1) double = 0
                titleStr(1,:) char   = ''
            end
            col = [0.2 0.5 0.9];
            figure('Name','Waveform Spectrum','Position',[100 100 800 400]);
            plot(f/1e6, Pdb, 'Color',col,'LineWidth',1.2); grid on;
            xlabel('Frequency [MHz]'); ylabel('PSD [dB, peak=0]');
            yline(-3, 'r--', '-3 dB', 'LineWidth',1);
            if B > 0
                xline([-B/2/1e6, B/2/1e6], 'k--', 'LineWidth',1);
            end
            if ~isempty(titleStr)
                title(titleStr);
            end
            ylim([-60 5]);
        end

    end

end
