classdef ReceiverChain < handle
% ReceiverChain  Analog/digital processing chain simulation.
%
%   Models the signal path:
%     RF → LNA → Mixer → IF BPF → VGA → AAF → ADC → DDC → Matched Filter
%
%   All filter design is windowed-sinc (no toolbox).  Filter orders are
%   derived automatically from the rate ratios and target attenuation.
%
%   CONSTRUCTION
%     ch = ReceiverChain(Name=Value, ...)
%
%   Optional Name-Value pairs
%     nBits      – ADC resolution; 0 = ideal    [bits]   (default: 0)
%     G_LNA_dB   – LNA power gain               [dB]     (default: 20)
%     fIF        – IF centre frequency          [Hz]     (default: 60e6)
%     fadc       – ADC sample rate              [Hz]     (default: 4*fIF)
%     fbb        – post-DDC baseband rate       [Hz]     (default: fadc/6)
%     Vfs        – ADC full-scale voltage       [V]      (default: 1)
%     headroom   – VGA RMS target as fraction   [-]      (default: 0.5)
%
%   METHODS
%     C = ch.process(s_rf, fan)           – run chain; no matched filter
%     C = ch.process(s_rf, fan, h_ref=h)  – run chain + pulse compression
%     ch.plot(C)                          – spectral + time-domain figures
%
%   CHAIN OUTPUT STRUCT  C
%     s_lna  – LNA output at fan  (real)
%     s_if   – IF BPF output at fan
%     s_vga  – VGA output at fan
%     s_aaf  – AAF output at fan
%     s_adc  – ADC output at fadc  (real, quantised)
%     s_bb   – DDC baseband at fbb  (complex)
%     y_mf   – matched filter output at fbb  ([] if no h_ref)
%     range  – range axis [m] corresponding to y_mf  ([] if no h_ref)
%     t_bb   – time axis [s] at fbb
%     fan, fadc, fbb, Gvga, d1, d2
%
%   EXAMPLE
%     B = 10e6;  Tp = 20e-6;  fIF = 60e6;  fan = 600e6;
%     t    = (0 : round(Tp*fan)-1).' / fan;
%     s_rf = cos(2*pi*fIF*t + pi*(B/Tp)*(t - Tp/2).^2);
%
%     ch   = ReceiverChain(nBits=12, fIF=fIF, fadc=150e6, fbb=25e6);
%     t_r  = (0 : round(Tp*ch.fbb)-1).' / ch.fbb;
%     h    = exp(1j*pi*(B/Tp)*(t_r - Tp/2).^2);
%     C    = ch.process(s_rf, fan, h_ref=h);
%     ch.plot(C);

    properties (SetAccess = private)
        nBits       % ADC resolution (0 = ideal)
        G_LNA_dB    % LNA power gain         [dB]
        fIF         % IF centre frequency    [Hz]
        fadc        % ADC sample rate        [Hz]
        fbb         % Baseband rate          [Hz]
        Vfs         % ADC full-scale voltage [V]
        headroom    % VGA RMS/Vfs target
    end

    properties (Dependent)
        dynamicRange_dB   % ADC SQNR: 6.02*nBits + 1.76 dB  (Inf if ideal)
    end

    properties (Constant, Access = private)
        c0 = 299792458;
    end

    % ==================================================================
    methods

        function obj = ReceiverChain(opts)
            arguments
                opts.nBits     (1,1) double {mustBeNonnegative} = 0
                opts.G_LNA_dB  (1,1) double                    = 20
                opts.fIF       (1,1) double {mustBePositive}   = 60e6
                opts.fadc      (1,1) double {mustBePositive}   = 0
                opts.fbb       (1,1) double {mustBePositive}   = 0
                opts.Vfs       (1,1) double {mustBePositive}   = 1
                opts.headroom  (1,1) double                    = 0.5
            end
            obj.nBits    = opts.nBits;
            obj.G_LNA_dB = opts.G_LNA_dB;
            obj.fIF      = opts.fIF;
            obj.Vfs      = opts.Vfs;
            obj.headroom = opts.headroom;
            obj.fadc     = ternary(opts.fadc > 0, opts.fadc, 4*opts.fIF);
            obj.fbb      = ternary(opts.fbb  > 0, opts.fbb,  obj.fadc/6);
        end

        function v = get.dynamicRange_dB(obj)
            if obj.nBits == 0, v = Inf; else, v = 6.02*obj.nBits + 1.76; end
        end

        % --------------------------------------------------------------
        function C = process(obj, s_rf, fan, opts)
        % process(s_rf, fan)            – chain simulation, no MF
        % process(s_rf, fan, h_ref=h)   – chain + FFT-based pulse compression
        %
        %   s_rf  : Nx1 real passband signal at IF, sample rate fan [Hz]
        %   h_ref : optional baseband MF reference at fbb (complex)
            arguments
                obj
                s_rf  (:,1) double
                fan   (1,1) double {mustBePositive}
                opts.NF_dB (1,1) double {mustBeNonnegative} = 0   % injected by RadarReceiver
                opts.T0    (1,1) double {mustBePositive}    = 290
                opts.h_ref (:,1) {mustBeNumeric}            = []
            end
            d1 = round(fan      / obj.fadc);
            d2 = round(obj.fadc / obj.fbb);
            assert(abs(d1 - fan/obj.fadc) < 1e-6, ...
                'ReceiverChain: fan must be an integer multiple of fadc.');
            assert(abs(d2 - obj.fadc/obj.fbb) < 1e-6, ...
                'ReceiverChain: fadc must be an integer multiple of fbb.');

            Nan  = numel(s_rf);
            kB   = 1.380649e-23;

            % ① LNA ───────────────────────────────────────────────────
            NF    = 10^(opts.NF_dB / 10);
            N0    = kB * opts.T0 * NF;
            Pn_an = N0 * fan;
            Gv    = 10^(obj.G_LNA_dB / 20);
            s_lna = Gv * (s_rf + sqrt(Pn_an)*randn(Nan,1));

            % ② Mixer (identity — LO translation done at synthesis) ───
            s_mix = s_lna;

            % ③ IF BPF ────────────────────────────────────────────────
            % Passband Bn centred at fIF; order for ~80 dB stopband.
            Bn   = opts.NF_dB;          % placeholder; real Bn comes from ReceiverNoise
            % Derive Bn from the signal itself if not injected properly:
            % Use ADC-derived passband as conservative bound: Bn_eff ~ fadc/d1
            Bn_eff = obj.fadc / d1;
            N_if  = Radar.ReceiverChain.filterOrder(fan, Bn_eff/2, 80);
            h_if  = Radar.ReceiverChain.bpf_wsinc(N_if, Bn_eff/2, obj.fIF, fan);
            s_if  = conv(s_mix, h_if, 'same');

            % ④ VGA ───────────────────────────────────────────────────
            Gvga  = obj.headroom * obj.Vfs / (sqrt(mean(s_if.^2)) + eps);
            s_vga = Gvga * s_if;

            % ⑤ AAF ───────────────────────────────────────────────────
            fc_aaf = 0.88 * obj.fadc / 2;
            N_aaf  = Radar.ReceiverChain.filterOrder(fan, fc_aaf, 80);
            h_aaf  = Radar.ReceiverChain.lpf_wsinc(N_aaf, fc_aaf, fan);
            s_aaf  = conv(s_vga, h_aaf, 'same');

            % ⑥ ADC ───────────────────────────────────────────────────
            s_down = s_aaf(1:d1:end);
            if obj.nBits > 0
                lsb   = 2*obj.Vfs / 2^obj.nBits;
                s_adc = max(min(round(s_down/lsb)*lsb, obj.Vfs-lsb), -obj.Vfs);
            else
                s_adc = s_down;
            end
            Nadc  = numel(s_adc);
            n_idx = (0:Nadc-1).';

            % ⑦ DDC ───────────────────────────────────────────────────
            s_demod = s_adc .* exp(-1j*2*pi*(obj.fIF/obj.fadc)*n_idx);
            fc_ddc  = 0.55 * obj.fbb / 2;
            N_ddc   = Radar.ReceiverChain.filterOrder(obj.fadc, fc_ddc, 60);
            h_ddc   = Radar.ReceiverChain.lpf_wsinc(N_ddc, fc_ddc, obj.fadc);
            s_lpf   = conv(s_demod, h_ddc, 'same');
            s_bb    = s_lpf(1:d2:end);
            Nbb     = numel(s_bb);
            t_bb    = (0:Nbb-1).' / obj.fbb;

            % ⑧ Matched filter ────────────────────────────────────────
            if ~isempty(opts.h_ref)
                Nfft  = 2^nextpow2(2*Nbb);
                Y     = fft(s_bb, Nfft) .* conj(fft(opts.h_ref(:), Nfft));
                y_mf  = ifft(Y);
                y_mf  = y_mf(1:Nbb);
                rng   = (0:Nbb-1).' / obj.fbb * (obj.c0/2);
            else
                y_mf  = [];
                rng   = [];
            end

            C = struct( ...
                's_lna', s_lna, 's_if', s_if, 's_vga', s_vga, ...
                's_aaf', s_aaf, 's_adc', s_adc, 's_bb', s_bb,  ...
                'y_mf',  y_mf,  'range', rng,  't_bb', t_bb,   ...
                'fan',   fan,   'fadc',  obj.fadc, 'fbb', obj.fbb, ...
                'Gvga',  Gvga,  'd1', d1, 'd2', d2);
        end

        % --------------------------------------------------------------
        function plot(obj, C)
        % plot(C)  — spectral evolution (Figure 1) + time snapshots (Figure 2).
            Nfig = 1024;

            % ── Figure 1: PSD at each stage ───────────────────────────
            figure('Name','ReceiverChain — Spectra','Position',[40 40 1320 860]);
            tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
            title(tl, sprintf( ...
                'Chain Spectra  |  f_{IF}=%.4g MHz  f_{adc}=%.4g MHz  f_{bb}=%.4g MHz', ...
                obj.fIF/1e6, obj.fadc/1e6, obj.fbb/1e6), 'FontSize',11);

            col = [0.2 0.5 0.9];

            nexttile;
            [f,P] = Radar.ReceiverChain.psd1(C.s_lna, C.fan, Nfig);
            plot(f/1e6, P, 'Color',col,'LineWidth',1); grid on;
            xlabel('Frequency [MHz]'); ylabel('PSD [dBW/Hz]');
            title('① LNA output');

            nexttile;
            [f,P] = Radar.ReceiverChain.psd1(C.s_if, C.fan, Nfig);
            plot(f/1e6, P, 'Color',col,'LineWidth',1); grid on;
            xlim([0 min(C.fan/2, 4*obj.fIF)/1e6]);
            xlabel('Frequency [MHz]'); ylabel('PSD [dBW/Hz]');
            title(sprintf('③ IF BPF  (centred %.4g MHz)', obj.fIF/1e6));

            nexttile;
            [f,P] = Radar.ReceiverChain.psd1(C.s_aaf, C.fan, Nfig);
            plot(f/1e6, P, 'Color',col,'LineWidth',1); grid on;
            xline(obj.fadc/2/1e6,'r--','f_s/2','LineWidth',1.2,...
                  'LabelOrientation','horizontal');
            xlim([0 min(C.fan/2, 1.5*obj.fadc)/1e6]);
            xlabel('Frequency [MHz]'); ylabel('PSD [dBW/Hz]');
            title(sprintf('⑤ AAF  (cutoff=%.4g MHz)', 0.88*obj.fadc/2/1e6));

            nexttile;
            [f,P] = Radar.ReceiverChain.psd1(real(C.s_adc), obj.fadc, Nfig);
            plot(f/1e6, P, 'Color',col,'LineWidth',1); grid on;
            xlim([0 obj.fadc/2/1e6]);
            xlabel('Frequency [MHz]'); ylabel('PSD [dBW/Hz]');
            if obj.nBits > 0
                title(sprintf('⑥ ADC  (%d-bit, SQNR≈%.0f dB)',obj.nBits,obj.dynamicRange_dB));
            else
                title('⑥ ADC  (ideal)');
            end

            nexttile;
            [f,P] = Radar.ReceiverChain.psd2(C.s_bb, obj.fbb, Nfig);
            plot(f/1e6, P, 'Color',col,'LineWidth',1); grid on;
            xlabel('Frequency [MHz]'); ylabel('PSD [dBW/Hz]');
            title(sprintf('⑦ DDC baseband  (f_{bb}=%.4g MHz)', obj.fbb/1e6));

            nexttile;
            if ~isempty(C.y_mf)
                y_db = 20*log10(abs(C.y_mf)/max(abs(C.y_mf))+eps);
                plot(C.range/1e3, y_db,'Color',col,'LineWidth',1);
                ylim([-60 5]); grid on;
                xlabel('Range [km]'); ylabel('Amplitude [dB]');
                title('⑧ Matched filter output');
            else
                axis off;
                text(0.5,0.5,'No matched filter reference supplied', ...
                    'Units','normalized','HorizontalAlignment','center');
            end

            % ── Figure 2: time snapshots ──────────────────────────────
            figure('Name','ReceiverChain — Time domain','Position',[80 80 1100 650]);
            tl2 = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
            title(tl2,'Time-domain snapshots','FontSize',11);

            Tshow = min(5e-6, numel(C.s_lna)/C.fan);
            t_an  = (0:numel(C.s_lna)-1).' / C.fan;
            n_ad  = (0:numel(C.s_adc)-1).';

            nexttile;
            mk = t_an <= Tshow;
            plot(t_an(mk)*1e6, C.s_lna(mk),'Color',col,'LineWidth',0.8);
            ylabel('Amplitude'); title('① LNA output'); grid on; xlabel('Time [µs]');

            nexttile;
            mk = n_ad/obj.fadc <= Tshow;
            stem(n_ad(mk)/obj.fadc*1e6, C.s_adc(mk), ...
                'Color',col,'MarkerSize',3,'LineWidth',0.8);
            yline([-obj.Vfs obj.Vfs],'r--','±V_{fs}');
            ylabel('Amplitude [V]');
            if obj.nBits > 0
                title(sprintf('⑥ ADC output (%d-bit)',obj.nBits));
            else
                title('⑥ ADC output (ideal)');
            end
            grid on; xlabel('Time [µs]');

            nexttile;
            mk = C.t_bb <= Tshow;
            if ~isempty(C.y_mf)
                plot(C.range(mk)/1e3, abs(C.y_mf(mk)),'Color',col,'LineWidth',1.2);
                ylabel('|y_{MF}|'); title('⑧ MF envelope'); xlabel('Range [km]');
            else
                plot(C.t_bb(mk)*1e6, abs(C.s_bb(mk)),'Color',col,'LineWidth',1);
                ylabel('|s_{bb}|'); title('⑦ Baseband envelope'); xlabel('Time [µs]');
            end
            grid on;
        end

        % --------------------------------------------------------------
        function disp(obj)
            fprintf('  ReceiverChain\n');
            fprintf('    G_LNA     = %.4g dB\n', obj.G_LNA_dB);
            fprintf('    fIF       = %s\n',      Radar.ReceiverChain.fmtHz(obj.fIF));
            fprintf('    fadc      = %s\n',      Radar.ReceiverChain.fmtHz(obj.fadc));
            fprintf('    fbb       = %s\n',      Radar.ReceiverChain.fmtHz(obj.fbb));
            fprintf('    Vfs       = %.4g V\n',  obj.Vfs);
            fprintf('    headroom  = %.4g\n',    obj.headroom);
            if obj.nBits > 0
                fprintf('    nBits     = %d  (SQNR≈%.0f dB)\n', ...
                    obj.nBits, obj.dynamicRange_dB);
            else
                fprintf('    nBits     = ideal\n');
            end
        end

    end   % public methods

    % ==================================================================
    methods (Static, Access = private)

        function N = filterOrder(fs, fc, atten_dB)
            df = min(fc, fs/2 - fc) / fs;
            N  = ceil(atten_dB / (22 * max(df, 1e-4)));
            N  = N + mod(N,2);
            N  = max(N, 64);
        end

        function h = lpf_wsinc(N, fc, fs)
            n  = (0:N).' - N/2;
            fn = fc / fs;
            w  = Radar.ReceiverChain.blackman(N+1);
            h  = 2*fn * Radar.ReceiverChain.sinc_fn(2*fn*n) .* w;
            h  = h / sum(h);
        end

        function h = bpf_wsinc(N, fc_half, f0, fs)
            n = (0:N).' - N/2;
            h = Radar.ReceiverChain.lpf_wsinc(N, fc_half, fs) .* ...
                2 .* cos(2*pi*(f0/fs)*n);
        end

        function w = blackman(N)
            n = (0:N-1).';
            w = 0.42 - 0.5*cos(2*pi*n/(N-1)) + 0.08*cos(4*pi*n/(N-1));
        end

        function y = sinc_fn(x)
            y     = ones(size(x));
            nz    = x ~= 0;
            y(nz) = sin(pi*x(nz)) ./ (pi*x(nz));
        end

        function [f, Pdb] = psd1(x, fs, N)
            X  = fft(x(:), N);
            P  = abs(X).^2 / (fs*N);
            P  = P(1:N/2+1);
            P(2:end-1) = 2*P(2:end-1);
            f  = (0:N/2).' * (fs/N);
            Pdb = 10*log10(P + eps);
        end

        function [f, Pdb] = psd2(x, fs, N)
            X   = fftshift(fft(x(:), N));
            P   = abs(X).^2 / (fs*N);
            f   = (-N/2:N/2-1).' * (fs/N);
            Pdb = 10*log10(P + eps);
        end

        function s = fmtHz(f)
            if     f >= 1e9, s = sprintf('%.4g GHz', f/1e9);
            elseif f >= 1e6, s = sprintf('%.4g MHz', f/1e6);
            elseif f >= 1e3, s = sprintf('%.4g kHz', f/1e3);
            else,            s = sprintf('%.4g Hz',  f);
            end
        end

    end

end

% ── File-scope helper ────────────────────────────────────────────────────
function v = ternary(cond, a, b)
    if cond, v = a; else, v = b; end
end
