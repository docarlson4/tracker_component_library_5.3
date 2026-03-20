% RadarReceiver_demo.m
% Feature-awareness script for RadarReceiver (and its RadarTransmitter
% coupling).  Run section-by-section or as a whole.  No toolbox required.
%
% Sections
%   1. Construction variants
%   2. Dependent properties & noise floor
%   3. noisePower — pre/post-filter bandwidth comparison
%   4. friisChain — LNA-dominated vs. lossy-front-end chains
%   5. snr() from a known received power
%   6. rangeSNR + maxRange via the radar range equation
%   7. Coherent vs. noncoherent integration gain & loss table
%   8. Albersheim curves — Pd vs. required SNR
%   9. Figures

close all; clc;
hdr = @(s) fprintf('\n%s\n%s\n', s, repmat('=', 1, numel(s)));

%% ========================================================================
hdr('1. CONSTRUCTION VARIANTS');
% =========================================================================

% Minimal (ideal ADC, no losses, T0 = 290 K)
rx0 = Radar.RadarReceiver(3, 1e6);
fprintf('rx0 (minimal):\n');  rx0

% Fully specified X-band receiver
rx = Radar.RadarReceiver(4, 5e6, losses_dB=1.5, T0=290, nBits=12);
fprintf('rx (X-band, 12-bit ADC):\n');  rx

% Cold front-end (space / cryogenic LNA scenario)
rxCold = Radar.RadarReceiver(0.5, 1e6, T0=50);
fprintf('rxCold (50 K front-end):\n');  rxCold

%% ========================================================================
hdr('2. DEPENDENT PROPERTIES & NOISE FLOOR');
% =========================================================================

fprintf('  Noise figure (linear F)       : %.4f\n',   rx.F);
fprintf('  Effective noise temp Te       : %.2f K\n', rx.Te);
fprintf('  kTBF at Bn = 5 MHz            : %.2f dBW\n', rx.kTBF_dBW);
fprintf('  MDS (incl. %.1f dB Rx loss)  : %.2f dBW\n', rx.losses_dB, rx.MDS_dBW);
fprintf('  ADC dynamic range (12-bit)    : %.2f dB\n',  rx.dynamicRange_dB);

% Show how MDS shifts with noise figure
NF_sweep = (1:0.5:10).';
MDS_sweep = arrayfun(@(nf) Radar.RadarReceiver(nf, 5e6, losses_dB=1.5).MDS_dBW, NF_sweep);
fprintf('\n  NF vs. MDS at Bn = 5 MHz, Lrx = 1.5 dB:\n');
fprintf('    NF = %4.1f dB  ->  MDS = %6.2f dBW\n', [NF_sweep, MDS_sweep].');

%% ========================================================================
hdr('3. noisePower — PRE/POST FILTER BANDWIDTH COMPARISON');
% =========================================================================

% Pre-detection (IF) bandwidth 50 MHz, matched-filter output ~1 MHz
BW_if   = 50e6;
BW_mf   =  1e6;
[Pn_if,   Pn_if_dBW]  = rx.noisePower(BW_if);
[Pn_mf,   Pn_mf_dBW]  = rx.noisePower(BW_mf);
fprintf('  Noise at IF  bandwidth (%3.0f MHz): %7.2f dBW  (%+.2f dB over MDS)\n', ...
    BW_if/1e6,  Pn_if_dBW,  Pn_if_dBW  - rx.MDS_dBW);
fprintf('  Noise at MF output  (%3.0f MHz): %7.2f dBW  (%+.2f dB over MDS)\n', ...
    BW_mf/1e6,  Pn_mf_dBW,  Pn_mf_dBW  - rx.MDS_dBW);
fprintf('  Noise reduction through MF filter : %.2f dB\n', Pn_if_dBW - Pn_mf_dBW);

%% ========================================================================
hdr('4. friisChain — CASCADE NOISE FIGURE');
% =========================================================================

% Chain A: LNA-first (well-designed)
%   LNA: NF=1.2 dB, G=25 dB
%   Cable: NF=0.8 dB (0.8 dB loss => G=-0.8 dB), passive
%   Mixer: NF=7 dB,  G=-6 dB
%   IF amp: NF=4 dB, G=30 dB
NF_A  = [1.2,  0.8,  7.0,  4.0];
G_A   = [25.0, -0.8, -6.0, 30.0];
NF_A_total = Radar.RadarReceiver.friisChain(NF_A, G_A);
fprintf('  Chain A (LNA-first): total NF = %.3f dB\n', NF_A_total);

% Chain B: lossy front-end (cable before LNA)
NF_B  = [0.8,  1.2,  7.0,  4.0];
G_B   = [-0.8, 25.0, -6.0, 30.0];
NF_B_total = Radar.RadarReceiver.friisChain(NF_B, G_B);
fprintf('  Chain B (cable-first): total NF = %.3f dB  (penalty = %.3f dB)\n', ...
    NF_B_total, NF_B_total - NF_A_total);

% Chain C: stage-by-stage cumulative NF (manual loop to show each stage)
fprintf('\n  Chain A cumulative NF per stage:\n');
for k = 1:numel(NF_A)
    fprintf('    After stage %d : %.4f dB\n', k, ...
        Radar.RadarReceiver.friisChain(NF_A(1:k), G_A(1:k)));
end

%% ========================================================================
hdr('5. snr() FROM KNOWN RECEIVED POWER');
% =========================================================================

% Suppose a coherent integration delivers Pr = -100 dBW at the port
Pr_vec_dBW = (-130:5:-90).';
Pr_vec_W   = 10.^(Pr_vec_dBW / 10);
[SNR_lin, SNR_dB_vec] = rx.snr(Pr_vec_W);
fprintf('  Received power -> SNR at Bn = %.0f MHz, NF = %.1f dB, Lrx = %.1f dB:\n', ...
    rx.Bn/1e6, rx.NF_dB, rx.losses_dB);
fprintf('    Pr = %5.0f dBW  ->  SNR = %6.2f dB\n', ...
    [Pr_vec_dBW, SNR_dB_vec].');

%% ========================================================================
hdr('6. rangeSNR + maxRange — RADAR RANGE EQUATION');
% =========================================================================

% X-band LFM system — Bn matched to post-compression bandwidth 1/tau
tx = Radar.RadarTransmitter(10e9, 100e3, 10e-6, 1e3, bandwidth=5e6, waveform='lfm');
rxLFM = Radar.RadarReceiver(4, 1/(10e-6), losses_dB=1.5);   % Bn = 1/tau = 100 kHz

fprintf('  Transmitter:\n');  tx
fprintf('  Receiver (Bn matched to tau):\n');  rxLFM

Gt_dBi   = 35;   Gr_dBi = 35;
sigma_1  = 1;    % 1 m^2 target

R_vec    = (10:10:500).' * 1e3;   % 10–500 km
[SNR_range, Pr_range] = rxLFM.rangeSNR(tx, Gt_dBi, Gr_dBi, sigma_1, R_vec);

Rmax_13  = rxLFM.maxRange(tx, Gt_dBi, Gr_dBi, sigma_1, 13);
Rmax_0   = rxLFM.maxRange(tx, Gt_dBi, Gr_dBi, sigma_1,  0);

fprintf('\n  Max range (SNRmin = 13 dB, sigma = 1 m^2): %.1f km\n', Rmax_13/1e3);
fprintf('  Max range (SNRmin =  0 dB, sigma = 1 m^2): %.1f km\n',  Rmax_0/1e3);

% RCS sweep: 0.01, 0.1, 1, 10, 100 m^2
sigma_vec = [0.01 0.1 1 10 100];
Rmax_sigma = arrayfun(@(s) rxLFM.maxRange(tx, Gt_dBi, Gr_dBi, s, 13), sigma_vec);
fprintf('\n  Rmax (SNRmin = 13 dB) vs. RCS:\n');
fprintf('    sigma = %6.2f m^2  ->  Rmax = %7.1f km\n', ...
    [sigma_vec; Rmax_sigma/1e3]);

%% ========================================================================
hdr('7. COHERENT / NONCOHERENT INTEGRATION GAIN & LOSS TABLE');
% =========================================================================

Pfa = 1e-6;   Pd = 0.9;
n_vec = [1 2 4 8 16 32 64 128 256];

fprintf('  Pfa = %.0e,  Pd = %.2f\n\n', Pfa, Pd);
fprintf('  %6s  %8s  %8s  %8s  %8s  %8s\n', ...
    'n', 'SNR1_dB', 'G_CI', 'G_NCI', 'L_int', 'SNRreq');
fprintf('  %s\n', repmat('-', 1, 60));

SNR1 = Radar.RadarReceiver.albersheim(Pfa, Pd, 1);
for n = n_vec
    G_CI  = rxLFM.coherentGain(n);
    G_NCI = rxLFM.noncoherentGain(Pfa, Pd, n);
    L_int = G_CI - G_NCI;
    SNRn  = Radar.RadarReceiver.albersheim(Pfa, Pd, n);
    fprintf('  %6d  %8.2f  %8.2f  %8.2f  %8.2f  %8.2f\n', ...
        n, SNR1, G_CI, G_NCI, L_int, SNRn);
end

%% ========================================================================
hdr('8. ALBERSHEIM CURVES — Pd vs. REQUIRED SNR');
% =========================================================================

Pfa_vals = [1e-3, 1e-6, 1e-9];
n_alb    = [1, 4, 16, 64];
Pd_axis  = (0.1:0.01:0.999).';

fprintf('  Single-pulse SNR1 required  [dB]  (Pd = 0.5, 0.9, 0.99):\n');
fprintf('  %8s', '');
fprintf('  Pd=%4.2f', [0.5, 0.9, 0.99]);
fprintf('\n');
for pf = Pfa_vals
    fprintf('  Pfa=%.0e', pf);
    for pd = [0.5, 0.9, 0.99]
        fprintf('  %7.2f', Radar.RadarReceiver.albersheim(pf, pd, 1));
    end
    fprintf('\n');
end

% Precompute Albersheim grid for plotting
SNR_alb = zeros(numel(Pd_axis), numel(n_alb));
for ki = 1:numel(n_alb)
    SNR_alb(:,ki) = arrayfun(@(p) Radar.RadarReceiver.albersheim(1e-6, p, n_alb(ki)), Pd_axis);
end

%% ========================================================================
hdr('9. FIGURES');
% =========================================================================

figure('Name','RadarReceiver — Feature Demo', ...
       'Position',[100 80 1000 600], 'Color','w');
tl = tiledlayout(2, 3, 'Padding','compact', 'TileSpacing','compact');

% --- Tile 1: MDS vs. NF --------------------------------------------------
nexttile;
NF_ax  = (0:0.25:12).';
MDS_ax = arrayfun(@(nf) Radar.RadarReceiver(nf, 5e6, losses_dB=1.5).MDS_dBW, NF_ax);
plot(NF_ax, MDS_ax, 'b-', LineWidth=1.5);
xline(rx.NF_dB,  '--r', sprintf('NF=%.0f dB', rx.NF_dB),  LabelVerticalAlignment='bottom');
yline(rx.MDS_dBW,'--r', sprintf('MDS=%.0f dBW',rx.MDS_dBW),LabelHorizontalAlignment='left');
xlabel('Noise Figure [dB]'); ylabel('MDS [dBW]');
title('Minimum Detectable Signal vs. NF');
grid on;

% --- Tile 2: Friis chain — per-stage cumulative NF -----------------------
nexttile;
stages = 1:numel(NF_A);
NF_cum_A = arrayfun(@(k) Radar.RadarReceiver.friisChain(NF_A(1:k), G_A(1:k)), stages);
NF_cum_B = arrayfun(@(k) Radar.RadarReceiver.friisChain(NF_B(1:k), G_B(1:k)), stages);
plot(stages, NF_cum_A, 'b-o', stages, NF_cum_B, 'r--s', LineWidth=1.5, MarkerSize=7);
legend('Chain A (LNA-first)', 'Chain B (cable-first)', Location='northwest');
xticks(stages);
xticklabels({'LNA/cable','LNA/cable','Mixer','IF amp'});
ylabel('Cumulative NF [dB]'); title('Friis Cascade — Cumulative NF'); grid on;

% --- Tile 3: SNR vs. received power --------------------------------------
nexttile;
Pr_fine  = (-140:1:-80).';
Pr_fine_W = 10.^(Pr_fine/10);
[~, SNR_fine] = rx.snr(Pr_fine_W);
plot(Pr_fine, SNR_fine, 'b-', LineWidth=1.5);
xline(rx.MDS_dBW, '--r', 'MDS', LabelVerticalAlignment='bottom');
yline(0, ':k');
xlabel('Received Power [dBW]'); ylabel('SNR [dB]');
title('SNR vs. Received Power  (Bn = 5 MHz)');
grid on;

% --- Tile 4: Range-SNR profile -------------------------------------------
nexttile;
plot(R_vec/1e3, SNR_range, 'b-', LineWidth=1.5);
yline(13, '--r', '13 dB threshold', LabelHorizontalAlignment='left');
yline(0,  ':k', '0 dB');
xline(Rmax_13/1e3, '--g', sprintf('%.0f km', Rmax_13/1e3), LabelVerticalAlignment='bottom');
xlabel('Range [km]'); ylabel('Single-Pulse SNR [dB]');
title('Range–SNR Profile  (\sigma = 1 m^2)');
xlim([0 500]); grid on;

% --- Tile 5: Rmax vs. RCS at various SNR thresholds ----------------------
nexttile;
sigma_ax  = logspace(-2, 3, 200).';
snr_thresholds = [0 10 13 20];
colors = lines(numel(snr_thresholds));
for ki = 1:numel(snr_thresholds)
    Rmax_ax = arrayfun(@(s) rxLFM.maxRange(tx, Gt_dBi, Gr_dBi, s, snr_thresholds(ki)), sigma_ax);
    loglog(sigma_ax, Rmax_ax/1e3, 'Color',colors(ki,:), LineWidth=1.5); hold on;
end
hold off;
legend(arrayfun(@(s) sprintf('SNR_{min} = %d dB', s), snr_thresholds, ...
       UniformOutput=false), Location='northwest');
xlabel('\sigma [m^2]'); ylabel('R_{max} [km]');
title('Maximum Range vs. RCS'); grid on;

% --- Tile 6: Albersheim — Pd vs. SNR (Pfa = 1e-6) -----------------------
nexttile;
colors6 = lines(numel(n_alb));
for ki = 1:numel(n_alb)
    semilogx(Pd_axis, SNR_alb(:,ki), 'Color',colors6(ki,:), LineWidth=1.5); hold on;
end
hold off;
legend(arrayfun(@(n) sprintf('n = %d', n), n_alb, UniformOutput=false), ...
       Location='southwest');
xlabel('P_d'); ylabel('Required SNR_1 [dB]');
title('Albersheim: P_d vs. SNR  (P_{fa} = 10^{-6})');
xlim([0.1 0.999]); grid on;

sgtitle('RadarReceiver — Feature Awareness', FontSize=14, FontWeight='bold');

fprintf('\n  All done.\n');
