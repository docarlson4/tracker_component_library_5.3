% RadarReceiver_demo.m
% Complete usage pattern for RadarReceiver.m
% Covers: construction, disp, link budget, Friis chain, Albersheim,
%         integration gains, range equation, and full processing chain.
%
% Waveform: LFM  B=10 MHz, Tp=20 µs (TBP=200), IF=60 MHz
% Rate hierarchy:
%   fan  = 600 MHz   analogue simulation (10× fIF, 60× B)
%   fadc = 150 MHz   ADC rate (fIF in 1st Nyquist zone, dec1=4)
%   fbb  =  25 MHz   post-DDC baseband (dec2=6, > 2B)

clear; close all; rng(42);

%% ── 1. Waveform / system parameters ─────────────────────────────────────
B       = 10e6;     % LFM bandwidth          [Hz]
Tp      = 20e-6;    % pulse width            [s]    → TBP = 200
fIF     = 60e6;     % IF centre frequency    [Hz]
fan     = 600e6;    % analogue sim rate      [Hz]
SNR_in  = 10;       % input SNR over B       [dB]   (pre-LNA)

%% ── 2. Construct RadarReceiver ───────────────────────────────────────────
rx = Radar.RadarReceiver( ...
    4,      ...     % NF_dB
    B,      ...     % Bn  (= IF filter BW = signal BW for matched receiver)
    losses_dB = 1.5,  ...
    nBits     = 12,   ...
    G_LNA_dB  = 20,   ...
    fIF       = fIF,  ...
    fadc      = 150e6, ...
    fbb       = 25e6,  ...
    Vfs       = 1.0,   ...
    headroom  = 0.5);

disp(rx);          % formatted property summary

%% ── 3. Friis noise chain ─────────────────────────────────────────────────
% LNA (NF=2 dB, G=20 dB) → Mixer (NF=8 dB, G=-6 dB) → IF amp (NF=5 dB, G=30 dB)
fprintf('\n── Friis noise chain ──\n');
NF_stages   = [2  8  5];
gain_stages = [20 -6 30];
NF_cas = Radar.RadarReceiver.friisChain(NF_stages, gain_stages);
fprintf('  Cascaded NF = %.4g dB\n', NF_cas);

%% ── 4. Albersheim: detection performance ────────────────────────────────
fprintf('\n── Albersheim (Pfa=1e-6, Pd=0.9) ──\n');
Pfa = 1e-6;  Pd = 0.9;

SNR_1    = Radar.RadarReceiver.albersheim(Pfa, Pd, 1);
SNR_16   = Radar.RadarReceiver.albersheim(Pfa, Pd, 16);
G_coh    = rx.coherentGain(16);
G_nci    = rx.noncoherentGain(Pfa, Pd, 16);
L_int    = G_coh - G_nci;

fprintf('  SNR req (n=1)  = %.4g dB\n', SNR_1);
fprintf('  SNR req (n=16) = %.4g dB\n', SNR_16);
fprintf('  Coherent gain  = %.4g dB\n', G_coh);
fprintf('  NCI gain       = %.4g dB\n', G_nci);
fprintf('  Integration loss = %.4g dB\n', L_int);

%% ── 5. Noise power and SNR at a given received power ────────────────────
fprintf('\n── Noise power ──\n');
[Pn_W, Pn_dBW] = rx.noisePower(B);
fprintf('  Noise power at Bn = %.4g dBW  (%.4g pW)\n', Pn_dBW, Pn_W*1e12);

Pr_test = 1e-13;    % −130 dBW received power example
[SNR_lin, SNR_dB_test] = rx.snr(Pr_test);
fprintf('  SNR at Pr=%.4g dBW : %.4g dB\n', 10*log10(Pr_test), SNR_dB_test);

%% ── 6. Range equation ────────────────────────────────────────────────────
% Requires a RadarTransmitter; build a minimal stand-in struct if the
% class is unavailable, or use the class directly.
%
% If RadarTransmitter is on path:
%   tx = Radar.RadarTransmitter(10e9, 100e3, Tp, 1e3, bandwidth=B, waveform='lfm');
%   R_vec  = (10:10:300).' * 1e3;
%   [SNR_range, Pr_range] = rx.rangeSNR(tx, 35, 35, 1, R_vec);
%   Rmax = rx.maxRange(tx, 35, 35, 1, SNR_1);
%   fprintf('\n── Range equation ──\n');
%   fprintf('  Max range (SNR=%.4g dB) = %.4g km\n', SNR_1, Rmax/1e3);

%% ── 7. Build analogue LFM at IF ─────────────────────────────────────────
fprintf('\n── Building analogue LFM signal at IF ──\n');
Nan = round(Tp * fan);
t   = (0:Nan-1).' / fan;

phi_tx  = pi * (B/Tp) * (t - Tp/2).^2;    % symmetric LFM phase
s_clean = cos(2*pi*fIF*t + phi_tx);        % unit-amplitude passband chirp

% Scale to desired pre-LNA SNR over B
kT   = 1.38e-23 * 290;
N0   = kT * 10^(rx.NF_dB/10);
Ps   = 10^(SNR_in/10) * N0 * B;
s_rf = sqrt(2*Ps) * s_clean;              % Nx1 real passband at IF

%% ── 8. Matched-filter reference at fbb ──────────────────────────────────
Nbb   = round(Tp * rx.fbb);
t_ref = (0:Nbb-1).' / rx.fbb;
h_ref = exp(1j * pi * (B/Tp) * (t_ref - Tp/2).^2);   % baseband complex chirp

%% ── 9. Run the processing chain ─────────────────────────────────────────
fprintf('── Running processing chain ──\n');
C = rx.processChain(s_rf, fan, h_ref=h_ref);

fprintf('  fan  → fadc  decimation : %d×\n', C.d1);
fprintf('  fadc → fbb   decimation : %d×\n', C.d2);
fprintf('  VGA gain                : %.4g dB\n', 20*log10(C.Gvga));
fprintf('  ADC RMS occupancy       : %.4g × Vfs\n', rms(C.s_adc)/rx.Vfs);
fprintf('  MF peak / range res     : %.4g m\n', rx.rangeRes);

%% ── 10. Plot the chain ───────────────────────────────────────────────────
rx.plotChain(C);

%% ── 11. Post-MF SNR estimate (sanity check) ──────────────────────────────
fprintf('\n── Post-MF SNR sanity check ──\n');
% Expected: SNR_in + compression gain (TBP dB) - NF - losses
TBP_dB      = 10*log10(B*Tp);
SNR_exp_dB  = SNR_in + TBP_dB - rx.NF_dB - rx.losses_dB;
fprintf('  Expected SNR (link budget) = %.4g dB\n', SNR_exp_dB);

% Measure from MF output: peak power vs noise floor in a guard region
y    = C.y_mf;
[pk, ipk] = max(abs(y));
guard = round(5 * rx.fbb / B);           % 5 range-resolution cells
noise_mask = true(size(y));
noise_mask(max(1,ipk-guard):min(end,ipk+guard)) = false;
SNR_meas_dB = 20*log10(pk) - 20*log10(rms(abs(y(noise_mask))) + eps);
fprintf('  Measured SNR  (MF output)  = %.4g dB\n', SNR_meas_dB);

%% ── 12. Spectrogram of baseband signal (pre-MF) ─────────────────────────
% Confirms linear FM sweep visible in the DDC output.
figure('Name','Baseband spectrogram','Position',[120 120 900 400]);
Nw   = 64;
hop  = 16;
win  = 0.5 - 0.5*cos(2*pi*(0:Nw-1)'/(Nw-1));   % Hann
[Sg, f_sg, t_sg] = local_stft(C.s_bb, rx.fbb, Nw, hop, win);
imagesc(t_sg*1e6, f_sg/1e6, 20*log10(abs(Sg)+eps));
axis xy; colorbar; colormap('jet');
clim([max(20*log10(abs(Sg(:))+eps))-60, max(20*log10(abs(Sg(:))+eps))]);
xlabel('Time [µs]'); ylabel('Frequency [MHz]');
title(sprintf('Baseband STFT spectrogram  (B=%g MHz, T_p=%g µs)', B/1e6, Tp*1e6));
yline([-B/2/1e6, B/2/1e6], 'w--', 'LineWidth', 1.2);

%% ── Local STFT ───────────────────────────────────────────────────────────
function [S, f, t] = local_stft(x, fs, Nw, hop, win)
% Minimal STFT, no toolbox.  Returns two-sided spectrum.
    x   = x(:);
    Nx  = numel(x);
    idx = 1 : hop : Nx - Nw + 1;
    S   = zeros(Nw, numel(idx));
    for k = 1:numel(idx)
        S(:,k) = fft(x(idx(k):idx(k)+Nw-1) .* win);
    end
    S = fftshift(S, 1);
    f = (-Nw/2:Nw/2-1).' * (fs/Nw);
    t = (idx + Nw/2 - 1).' / fs;
end
