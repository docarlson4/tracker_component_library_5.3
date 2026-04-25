% RadarSystem_demo.m
% Unified usage pattern for the refactored transmitter + receiver hierarchy.
%
%   TransmitterParams / TransmitterWaveform / RadarTransmitter
%   ReceiverNoise    / ReceiverChain       / ReceiverDetection / RadarReceiver
%
% Key design: signal construction is explicit and external to both classes.
% The waveform builders (TransmitterWaveform static methods) and the chain
% (ReceiverChain.process) are wired together in the script, not hidden inside
% monolithic class methods.
clc
clear; close all; rng(42);
dbstop if error

%% ═══════════════════════════════════════════════════════════════════════
%% 1.  SYSTEM PARAMETERS
%% ═══════════════════════════════════════════════════════════════════════
fc   = 10e9;      % carrier frequency     [Hz]   X-band
B    = 10e6;      % LFM bandwidth         [Hz]
Tp   = 20e-6;     % pulse width           [s]    TBP = 200
PRF  = 1e3;       % pulse repetition freq [Hz]
Pt   = 100e3;     % peak power            [W]
fan  = 600e6;     % analogue sim rate     [Hz]   10× fIF, 60× B
fIF  = 60e6;      % IF centre freq        [Hz]
SNR_in = 10;      % pre-LNA input SNR     [dB]

%% ═══════════════════════════════════════════════════════════════════════
%% 2.  BUILD SUB-OBJECTS DIRECTLY
%% ═══════════════════════════════════════════════════════════════════════

% ── Transmitter ───────────────────────────────────────────────────────
txp = Radar.TransmitterParams(fc, Pt, Tp, PRF, bandwidth=B, waveform='lfm', losses_dB=0.5);
disp(txp);

% ── Receiver ──────────────────────────────────────────────────────────
noise = Radar.ReceiverNoise(4, B, losses_dB=1.5);
chain = Radar.ReceiverChain(nBits=12, G_LNA_dB=20, fIF=fIF, fadc=150e6, fbb=25e6);
disp(noise);
disp(chain);

%% ═══════════════════════════════════════════════════════════════════════
%% 3.  COMPOSE FACADES
%% ═══════════════════════════════════════════════════════════════════════
tx = Radar.RadarTransmitter(txp);         % pass pre-built params
rx = Radar.RadarReceiver(noise, chain);   % pass pre-built sub-objects
disp(tx);
disp(rx);

%% ═══════════════════════════════════════════════════════════════════════
%% 4.  TRANSMITTER: WAVEFORM CONSTRUCTION (explicit, outside the class)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── External waveform construction ──\n');

% ① Time axis for the analogue signal
Nan = round(Tp * fan);
t   = (0:Nan-1).' / fan;

% ② Build LFM envelope directly — no transmitter object required
s_base = Radar.TransmitterWaveform.lfm(t, Tp, B);    % unit amplitude, baseband

% ③ Scale to desired pre-LNA SNR over B
kT   = 1.38e-23 * 290;
N0   = kT * 10^(noise.NF_dB/10);
Ps   = 10^(SNR_in/10) * N0 * B;
s_if = sqrt(2*Ps) * real(s_base .* exp(1j*2*pi*fIF*t));   % passband at IF

fprintf('  Signal RMS  = %.4g V\n', rms(s_if));
fprintf('  Signal power = %.4g dBW\n', 10*log10(mean(s_if.^2)));

% ④ Alternatively: use the facade (dispatches to TransmitterWaveform.pulse)
%    s_tx is the complex envelope; translate to passband manually for chain
[s_env, t_tx] = tx.generatePulse(fan);
s_if_alt = sqrt(2*Ps) / rms(real(s_env)) * real(s_env .* exp(1j*2*pi*fIF*t_tx));
fprintf('  Both methods agree: %d\n', norm(s_if - s_if_alt) < 1e-10*norm(s_if));

%% ═══════════════════════════════════════════════════════════════════════
%% 5.  TRANSMITTER: SPECTRUM
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Waveform spectrum ──\n');
fs_spec = 100e6;
[Pdb, f_sp, BW_3dB] = Radar.TransmitterWaveform.spectrum(s_base, fs_spec);
fprintf('  3-dB bandwidth = %.4g MHz  (expected ≈ %.4g MHz)\n', BW_3dB/1e6, B/1e6);
Radar.TransmitterWaveform.plotSpectrum(Pdb(:), f_sp(:), B, ...
    sprintf('LFM spectrum  B=%g MHz, T_p=%g µs', B/1e6, Tp*1e6));

%% ═══════════════════════════════════════════════════════════════════════
%% 6.  TRANSMITTER: AMBIGUITY FUNCTION
%% ═══════════════════════════════════════════════════════════════════════
fprintf('── Ambiguity function ──\n');
fs_af = 50e6;
t_af  = (0:1/fs_af:Tp-1/fs_af).';
s_af  = Radar.TransmitterWaveform.lfm(t_af, Tp, B);                  % direct builder
[chi, tau_ax, fd_ax] = Radar.TransmitterWaveform.ambiguityFunction(s_af, fs_af, 256);
fprintf('  Delay axis: %.4g µs to %.4g µs\n', tau_ax(1)*1e6, tau_ax(end)*1e6);
Radar.TransmitterWaveform.plotAmbiguity(chi, tau_ax, fd_ax, ...
    sprintf('LFM Ambiguity  B=%g MHz, T_p=%g µs (TBP=%g)', B/1e6, Tp*1e6, B*Tp));

%% ═══════════════════════════════════════════════════════════════════════
%% 7.  RECEIVER: DETECTION THRESHOLDS
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── ReceiverDetection ──\n');
Pfa = 1e-6;  Pd = 0.9;  nInt = 16;

SNR_1  = Radar.ReceiverDetection.albersheim(Pfa, Pd, 1);
SNR_n  = Radar.ReceiverDetection.albersheim(Pfa, Pd, nInt);
G_ci   = Radar.ReceiverDetection.coherentGain(nInt);
G_nci  = Radar.ReceiverDetection.noncoherentGain(Pfa, Pd, nInt);
L_int  = Radar.ReceiverDetection.integrationLoss(Pfa, Pd, nInt);

fprintf('  Required SNR (n=1)   = %.4g dB\n', SNR_1);
fprintf('  Required SNR (n=%d)  = %.4g dB\n', nInt, SNR_n);
fprintf('  Coherent gain        = %.4g dB\n', G_ci);
fprintf('  NCI gain             = %.4g dB\n', G_nci);
fprintf('  Integration loss     = %.4g dB\n', L_int);

%% ═══════════════════════════════════════════════════════════════════════
%% 8.  RANGE EQUATION
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Range equation ──\n');
Gt_dBi = 35;   Gr_dBi = 35;   sigma = 1;   % [dBi dBi m^2]

R_vec  = (10:10:300).' * 1e3;
[SNR_r, Pr_r] = rx.rangeSNR(tx, Gt_dBi, Gr_dBi, sigma, R_vec);
Rmax   = rx.maxRange(tx, Gt_dBi, Gr_dBi, sigma, SNR_1);

fprintf('  Max range (Pd=%.2f, Pfa=%.0e) = %.4g km\n', Pd, Pfa, Rmax/1e3);

figure('Name','Range equation','Position',[50 50 900 400]);
tl_r = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl_r,'Range equation','FontSize',11);
col = [0.2 0.5 0.9];

nexttile;
plot(R_vec/1e3, SNR_r,'Color',col,'LineWidth',1.5); grid on;
yline(SNR_1,'r--',sprintf('%.0f dB (Albersheim)',SNR_1),'LineWidth',1.2);
xline(Rmax/1e3,'k--',sprintf('R_{max}=%.0f km',Rmax/1e3),'LineWidth',1);
xlabel('Range [km]'); ylabel('SNR [dB]'); title('Single-pulse SNR vs Range');

nexttile;
plot(R_vec/1e3, Pr_r,'Color',col,'LineWidth',1.5); grid on;
yline(noise.MDS_dBW,'r--','MDS','LineWidth',1.2);
xlabel('Range [km]'); ylabel('P_r [dBW]'); title('Received power vs Range');

%% ═══════════════════════════════════════════════════════════════════════
%% 9.  MF REFERENCE — built explicitly at fbb
%% ═══════════════════════════════════════════════════════════════════════
Nbb   = round(Tp * chain.fbb);
t_ref = (0:Nbb-1).' / chain.fbb;
h_ref = Radar.TransmitterWaveform.lfm(t_ref, Tp, B);   % baseband complex chirp at fbb

%% ═══════════════════════════════════════════════════════════════════════
%% 10. RECEIVER PROCESSING CHAIN
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Running receiver chain ──\n');

% Via facade: injects NF_dB and T0 automatically
C = rx.process(s_if, fan, h_ref=h_ref);

% Alternatively — chain directly, NF_dB and T0 explicit:
% C = chain.process(s_if, fan, NF_dB=noise.NF_dB, T0=noise.T0, h_ref=h_ref);

fprintf('  dec1 (fan → fadc) = %d×\n', C.d1);
fprintf('  dec2 (fadc → fbb) = %d×\n', C.d2);
fprintf('  VGA gain          = %.4g dB\n', 20*log10(C.Gvga));
fprintf('  ADC occupancy     = %.4g × V_fs\n', rms(C.s_adc)/chain.Vfs);

rx.plot(C);   % delegates to chain.plot(C)

%% ═══════════════════════════════════════════════════════════════════════
%% 11. POST-MF SNR SANITY CHECK
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Post-MF SNR sanity check ──\n');
TBP_dB     = 10*log10(B*Tp);
SNR_exp_dB = SNR_in + TBP_dB - noise.NF_dB - noise.losses_dB;
fprintf('  Expected (link budget) = %.4g dB\n', SNR_exp_dB);

y = C.y_mf;
[pk, ipk] = max(abs(y));
guard      = round(5 * chain.fbb / B);
mask       = true(size(y));
mask(max(1,ipk-guard):min(end,ipk+guard)) = false;
SNR_meas   = 20*log10(pk) - 20*log10(rms(abs(y(mask)))+eps);
fprintf('  Measured (MF output)   = %.4g dB\n', SNR_meas);

%% ═══════════════════════════════════════════════════════════════════════
%% 12. DDC BASEBAND SPECTROGRAM
%% ═══════════════════════════════════════════════════════════════════════
figure('Name','Baseband spectrogram','Position',[120 120 900 400]);
Nw  = 64;  hop = 16;
win = 0.5 - 0.5*cos(2*pi*(0:Nw-1)'/(Nw-1));
[Sg, f_sg, t_sg] = local_stft(C.s_bb, chain.fbb, Nw, hop, win);
imagesc(t_sg*1e6, f_sg/1e6, 20*log10(abs(Sg)+eps));
axis xy; colorbar; colormap('jet');
clim([max(20*log10(abs(Sg(:))+eps))-60, max(20*log10(abs(Sg(:))+eps))]);
xlabel('Time [µs]'); ylabel('Frequency [MHz]');
title(sprintf('DDC baseband STFT  (B=%g MHz, T_p=%g µs)', B/1e6, Tp*1e6));
yline([-B/2/1e6, B/2/1e6],'w--','LineWidth',1.2);

%% ═══════════════════════════════════════════════════════════════════════
%% 13. BARKER-13 PHASE CODE  (using builders directly, no facade)
%% ═══════════════════════════════════════════════════════════════════════
fprintf('\n── Barker-13 waveform (direct builder) ──\n');
code_b13 = pi * [0 0 0 0 0 1 1 0 0 1 0 1 0].';
fs_b     = 13e6;                                 % 1 sample per chip
N_b      = 13;
s_b13    = Radar.TransmitterWaveform.phaseCode(N_b, code_b13);
h_b13    = Radar.TransmitterWaveform.matchedFilter(s_b13);
y_b13    = conv(s_b13, h_b13);

figure('Name','Barker-13','Position',[140 140 900 350]);
tl_b = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
title(tl_b,'Barker-13 phase code','FontSize',11);

nexttile;
stem(0:N_b-1, angle(s_b13)/pi,'Color',col,'LineWidth',1.5,'MarkerSize',5);
xlabel('Chip index'); ylabel('Phase [π rad]'); grid on;
title('Phase code');

nexttile;
n_corr = -(N_b-1):N_b-1;
plot(n_corr, abs(y_b13),'Color',col,'LineWidth',1.5); grid on;
xlabel('Lag [chips]'); ylabel('|Autocorrelation|');
title(sprintf('Compression output  (PSLR = %.4g dB)', ...
    20*log10(1/13)));   % Barker-13 PSLR = -22.3 dB

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
