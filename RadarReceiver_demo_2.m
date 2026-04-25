% RadarReceiver_demo.m
% Usage pattern for the refactored receiver class hierarchy:
%   ReceiverNoise / ReceiverChain / ReceiverDetection / RadarReceiver

clear; close all; rng(42);

%% ── 1. Waveform / system parameters ─────────────────────────────────────
B    = 10e6;      % LFM bandwidth      [Hz]
Tp   = 20e-6;     % pulse width        [s]   → TBP = 200
fIF  = 60e6;      % IF centre          [Hz]
fan  = 600e6;     % analogue sim rate  [Hz]

%% ── 2. Build sub-objects independently ──────────────────────────────────

Bn = 2*B;
noise = Radar.ReceiverNoise(4, Bn, losses_dB=1.5, T0=290);
chain = Radar.ReceiverChain(nBits=12, G_LNA_dB=20, ...
                      fIF=fIF, fadc=150e6, fbb=25e6, ...
                      Vfs=1.0, headroom=0.5);

disp(noise);
disp(chain);

%% ── 3. Compose into facade ───────────────────────────────────────────────
rx = Radar.RadarReceiver(noise, chain);
disp(rx);

%% ── 4. Noise chain ───────────────────────────────────────────────────────
fprintf('\n── ReceiverNoise ──\n');

% Friis: LNA (NF=2 dB, G=20) → Mixer (8 dB, -6 dB) → IF amp (5 dB, 30 dB)
NF_cas = Radar.ReceiverNoise.friisChain([2 8 5], [20 -6 30]);
fprintf('  Cascaded NF   = %.4g dB\n', NF_cas);

[Pn_W, Pn_dBW] = noise.noisePower(Bn);
fprintf('  Noise power   = %.4g dBW  (%.4g pW)\n', Pn_dBW, Pn_W*1e12);

Pr_test = 1e-13;
[~, SNR_test] = noise.snr(Pr_test);
fprintf('  SNR at Pr=%.0f dBW : %.4g dB\n', 10*log10(Pr_test), SNR_test);

%% ── 5. Detection ─────────────────────────────────────────────────────────
fprintf('\n── ReceiverDetection ──\n');
Pfa = 1e-6;  Pd = 0.9;  n = 16;

SNR_1  = Radar.ReceiverDetection.albersheim(Pfa, Pd, 1);
SNR_n  = Radar.ReceiverDetection.albersheim(Pfa, Pd, n);
G_ci   = Radar.ReceiverDetection.coherentGain(n);
G_nci  = Radar.ReceiverDetection.noncoherentGain(Pfa, Pd, n);
L_int  = Radar.ReceiverDetection.integrationLoss(Pfa, Pd, n);

fprintf('  SNR req (n=1)     = %.4g dB\n', SNR_1);
fprintf('  SNR req (n=%d)    = %.4g dB\n', n, SNR_n);
fprintf('  Coherent gain     = %.4g dB\n', G_ci);
fprintf('  NCI gain          = %.4g dB\n', G_nci);
fprintf('  Integration loss  = %.4g dB\n', L_int);

% Also accessible through facade:
assert(rx.detect.albersheim(Pfa, Pd, 1) == SNR_1);

%% ── 6. Range equation (requires RadarTransmitter on path) ───────────────
kT     = 1.38e-23 * 290;
N0     = kT * 10^(noise.NF_dB/10);
SNR_in = 10;  % [dB] input SNR before LNA
Ps     = 10^(SNR_in/10) * N0 * B;

tx     = Radar.RadarTransmitter(fIF, Ps, Tp, 1e3, bandwidth=B, waveform='lfm');
R_vec  = (10:10:300).' * 1e3;
[SNR_range, ~] = rx.rangeSNR(tx, 35, 35, 1, R_vec);
Rmax = rx.maxRange(tx, 35, 35, 1, SNR_1);
fprintf('\n  Max range (SNR=%.4g dB) = %.4g km\n', SNR_1, Rmax/1e3);

[s, t] = tx.generatePulse(fan);

% ── 7. Build analogue LFM signal ─────────────────────────────────────────
fprintf('\n── Building LFM signal ──\n');
Nan = round(Tp * fan);
t   = (0:Nan-1).' / fan;

s_rf   = sqrt(2*Ps) * cos(2*pi*fIF*t + pi*(B/Tp)*(t - Tp/2).^2);
figure,plot(t,sqrt(2*Ps) * real(s),t,s_rf)

%% ── 8. MF reference at baseband / fbb ───────────────────────────────────
Nbb   = round(Tp * chain.fbb);
t_ref = (0:Nbb-1).' / chain.fbb;
h_ref = exp(1j * pi * (B/Tp) * (t_ref - Tp/2).^2);

%% ── 9a. Run chain via facade (NF_dB and T0 injected automatically) ───────
fprintf('── Running chain via RadarReceiver facade ──\n');
C = rx.process(s_rf, fan, h_ref=h_ref);
fprintf('  dec1 (fan→fadc) = %d×\n', C.d1);
fprintf('  dec2 (fadc→fbb) = %d×\n', C.d2);
fprintf('  VGA gain        = %.4g dB\n', 20*log10(C.Gvga));

%% ── 9b. OR run chain directly (must pass NF_dB, T0 explicitly) ──────────
% C2 = chain.process(s_rf, fan, NF_dB=noise.NF_dB, T0=noise.T0, h_ref=h_ref);

%% ── 10. Plot ─────────────────────────────────────────────────────────────
rx.plot(C);      % delegates to chain.plot(C)
% chain.plot(C); % equivalent direct call

%% ── 11. Post-MF SNR sanity check ────────────────────────────────────────
fprintf('\n── Post-MF SNR sanity check ──\n');
TBP_dB     = 10*log10(B*Tp);
SNR_exp_dB = SNR_in + TBP_dB - noise.NF_dB - noise.losses_dB;
fprintf('  Expected (budget) = %.4g dB\n', SNR_exp_dB);

y = C.y_mf;
[pk, ipk] = max(abs(y));
guard      = round(5 * chain.fbb / B);
mask       = true(size(y));
mask(max(1,ipk-guard):min(end,ipk+guard)) = false;
SNR_meas   = 20*log10(pk) - 20*log10(rms(abs(y(mask)))+eps);
fprintf('  Measured (MF out) = %.4g dB\n', SNR_meas);

%% ── 12. Baseband spectrogram ─────────────────────────────────────────────
figure('Name','Baseband spectrogram','Position',[120 120 900 400]);
Nw  = 64;  hop = 16;
win = 0.5 - 0.5*cos(2*pi*(0:Nw-1)'/(Nw-1));
[Sg, f_sg, t_sg] = local_stft(C.s_bb, chain.fbb, Nw, hop, win);
imagesc(t_sg*1e6, f_sg/1e6, 20*log10(abs(Sg)+eps));
axis xy; colorbar; colormap('jet');
clim([max(20*log10(abs(Sg(:))+eps))-60, max(20*log10(abs(Sg(:))+eps))]);
xlabel('Time [µs]'); ylabel('Frequency [MHz]');
title(sprintf('DDC baseband STFT  (B=%g MHz, T_p=%g µs)', B/1e6, Tp*1e6));
yline([-B/2/1e6, B/2/1e6], 'w--', 'LineWidth',1.2);

%% ── Local STFT (no toolbox) ──────────────────────────────────────────────
function [S, f, t] = local_stft(x, fs, Nw, hop, win)
    x   = x(:);
    idx = 1:hop:numel(x)-Nw+1;
    S   = zeros(Nw, numel(idx));
    for k = 1:numel(idx)
        S(:,k) = fft(x(idx(k):idx(k)+Nw-1) .* win);
    end
    S = fftshift(S,1);
    f = (-Nw/2:Nw/2-1).' * (fs/Nw);
    t = (idx + Nw/2 - 1).' / fs;
end
