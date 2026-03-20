classdef RadarReceiver < handle
% RadarReceiver  Parametric radar receiver: noise, sensitivity, range
%               equation, and pulse-integration analysis.
%
%   No toolbox dependencies.  Requires R2021a+.
%
%   CONSTRUCTION
%     rx = RadarReceiver(NF_dB, Bn)
%     rx = RadarReceiver(NF_dB, Bn, Name=Value, ...)
%
%   Required arguments
%     NF_dB  – receiver noise figure           [dB]
%     Bn     – noise (pre-detection) bandwidth [Hz]
%
%   Optional Name-Value pairs
%     losses_dB  – ohmic receive-chain losses  [dB]    (default: 0)
%     T0         – reference temperature       [K]     (default: 290)
%     nBits      – ADC resolution; 0 = ideal   [bits]  (default: 0)
%
%   EXAMPLES
%     % Sensitivity and noise chain
%     rx = RadarReceiver(4, 5e6, losses_dB=1.5, nBits=12);
%     rx
%
%     % Friis chain: LNA (2 dB NF, 20 dB gain) -> mixer (8 dB NF, -6 dB) -> IF amp (5 dB NF, 30 dB)
%     NF_total = RadarReceiver.friisChain([2 8 5], [20 -6 30])
%
%     % Range equation with a RadarTransmitter
%     tx     = RadarTransmitter(10e9, 100e3, 10e-6, 1e3, bandwidth=5e6, waveform='lfm');
%     rx     = RadarReceiver(4, 5e6);
%     R_vec  = (10:10:200).' * 1e3;
%     SNR_dB = rx.rangeSNR(tx, 35, 35, 1, R_vec);    % Gt=Gr=35 dBi, RCS=1 m^2
%     Rmax   = rx.maxRange(tx, 35, 35, 1, 13)         % 13 dB SNR threshold [m]
%
%     % Albersheim: required per-pulse SNR, Pd=0.9, Pfa=1e-6, 16-pulse NCI
%     SNR_req  = RadarReceiver.albersheim(1e-6, 0.9, 16)    % [dB]
%     G_NCI    = rx.noncoherentGain(1e-6, 0.9, 16)          % [dB]
%     L_int    = rx.coherentGain(16) - G_NCI                % integration loss [dB]

    % ------------------------------------------------------------------
    properties (SetAccess = private)
        NF_dB       % Receiver noise figure [dB]
        Bn          % Noise bandwidth [Hz]
        losses_dB   % Receive-chain ohmic losses [dB]
        T0          % Reference noise temperature [K]
        nBits       % ADC resolution (0 = ideal)
    end

    % ------------------------------------------------------------------
    properties (Dependent)
        F               % Noise figure (linear)
        Te              % Effective noise temperature  Te = T0*(F-1)  [K]
        kTBF_W          % Thermal noise power at Bn  [W]
        kTBF_dBW        % Thermal noise power at Bn  [dBW]
        MDS_dBW         % Minimum detectable signal (SNR=0 dB, incl. Rx losses) [dBW]
        dynamicRange_dB % ADC SFDR: 6.02*nBits + 1.76 dB  (Inf for ideal)
    end

    % ------------------------------------------------------------------
    properties (Constant, Access = private)
        kB = 1.380649e-23;   % Boltzmann constant [J/K]
        c0 = 299792458;      % Speed of light [m/s]
    end

    % ==================================================================
    methods

        % --------------------------------------------------------------
        function obj = RadarReceiver(NF_dB, Bn, options)
            arguments
                NF_dB      (1,1) double {mustBeNonnegative}
                Bn         (1,1) double {mustBePositive}
                options.losses_dB (1,1) double {mustBeNonnegative} = 0
                options.T0        (1,1) double {mustBePositive}    = 290
                options.nBits     (1,1) double {mustBeNonnegative, mustBeInteger} = 0
            end
            obj.NF_dB     = NF_dB;
            obj.Bn        = Bn;
            obj.losses_dB = options.losses_dB;
            obj.T0        = options.T0;
            obj.nBits     = options.nBits;
        end

        % --------------------------------------------------------------
        % Dependent getters
        % --------------------------------------------------------------
        function v = get.F(obj)
            v = 10^(obj.NF_dB / 10);
        end
        function v = get.Te(obj)
            v = obj.T0 * (obj.F - 1);
        end
        function v = get.kTBF_W(obj)
            v = obj.kB * obj.T0 * obj.Bn * obj.F;
        end
        function v = get.kTBF_dBW(obj)
            v = 10*log10(obj.kTBF_W);
        end
        function v = get.MDS_dBW(obj)
            v = obj.kTBF_dBW + obj.losses_dB;
        end
        function v = get.dynamicRange_dB(obj)
            if obj.nBits == 0
                v = Inf;
            else
                v = 6.02 * obj.nBits + 1.76;
            end
        end

        % --------------------------------------------------------------
        function [Pn_W, Pn_dBW] = noisePower(obj, BW)
        % noisePower(BW)  ->  thermal noise power [W] and [dBW] at bandwidth BW.
        %
        %   Useful for comparing noise in pre- and post-filter bandwidths,
        %   or for computing noise figure of sub-bands.
            arguments
                obj
                BW (1,1) double {mustBePositive}
            end
            Pn_W   = obj.kB * obj.T0 * BW * obj.F;
            Pn_dBW = 10*log10(Pn_W);
        end

        % --------------------------------------------------------------
        function [SNR, SNR_dB] = snr(obj, Pr_W)
        % snr(Pr_W)  ->  SNR [linear] and [dB] at noise bandwidth Bn.
        %
        %   Pr_W may be a vector (e.g. swept over range or time).
        %   Includes receive-chain losses.
            arguments
                obj
                Pr_W (:,1) double {mustBePositive}
            end
            Lrx    = 10^(obj.losses_dB / 10);
            SNR    = Pr_W ./ (obj.kTBF_W * Lrx);
            SNR_dB = 10*log10(SNR);
        end

        % --------------------------------------------------------------
        function [SNR_dB, Pr_dBW] = rangeSNR(obj, tx, Gt_dBi, Gr_dBi, sigma_m2, R_m)
        % rangeSNR(tx, Gt_dBi, Gr_dBi, sigma_m2, R_m)
        %
        %   Monostatic radar range equation evaluated at range(s) R_m.
        %
        %   tx       : RadarTransmitter object (supplies Pt, lambda, losses_dB)
        %   Gt_dBi   : transmit antenna gain       [dBi]
        %   Gr_dBi   : receive  antenna gain       [dBi]
        %   sigma_m2 : target RCS                  [m^2]  (scalar)
        %   R_m      : slant range                 [m]    (Nx1 vector or scalar)
        %
        %   System loss = tx.losses_dB + rx.losses_dB (combined one-way
        %   ohmic / distribution losses; do not double-count).
        %
        %   Returns
        %     SNR_dB  : Nx1  single-pulse SNR at receiver output   [dB]
        %     Pr_dBW  : Nx1  received signal power at antenna port [dBW]
            arguments
                obj
                tx           Radar.RadarTransmitter
                Gt_dBi  (1,1) double
                Gr_dBi  (1,1) double
                sigma_m2(1,1) double {mustBePositive}
                R_m     (:,1) double {mustBePositive}
            end
            Gt   = 10^(Gt_dBi  / 10);
            Gr   = 10^(Gr_dBi  / 10);
            Ltot = 10^((tx.losses_dB + obj.losses_dB) / 10);
            lam  = tx.lambda;

            num    = tx.Pt * Gt * Gr * lam^2 * sigma_m2;
            denom  = (4*pi)^3 * R_m.^4 * Ltot;

            Pr_W   = num ./ denom;
            Pr_dBW = 10*log10(Pr_W);
            SNR_dB = Pr_dBW - obj.kTBF_dBW;
        end

        % --------------------------------------------------------------
        function R_m = maxRange(obj, tx, Gt_dBi, Gr_dBi, sigma_m2, SNRmin_dB)
        % maxRange(tx, Gt_dBi, Gr_dBi, sigma_m2, SNRmin_dB)  ->  [m]
        %
        %   Closed-form radar range equation solved for detection range.
        %
        %   SNRmin_dB is the required single-pulse output SNR.  For a
        %   pulse-compressed waveform, Bn should equal the post-compression
        %   noise bandwidth (≈ 1/tau for time-bandwidth-matched filter).
            arguments
                obj
                tx              Radar.RadarTransmitter
                Gt_dBi    (1,1) double
                Gr_dBi    (1,1) double
                sigma_m2  (1,1) double {mustBePositive}
                SNRmin_dB (1,1) double
            end
            Gt     = 10^(Gt_dBi    / 10);
            Gr     = 10^(Gr_dBi    / 10);
            SNRmin = 10^(SNRmin_dB / 10);
            Ltot   = 10^((tx.losses_dB + obj.losses_dB) / 10);

            num  = tx.Pt * Gt * Gr * tx.lambda^2 * sigma_m2;
            den  = (4*pi)^3 * obj.kTBF_W * Ltot * SNRmin;
            R_m  = (num / den)^0.25;
        end

        % --------------------------------------------------------------
        function G_dB = coherentGain(~, nPulses)
        % coherentGain(nPulses)  ->  coherent integration gain [dB].
        %
        %   Exact: G_CI = 10*log10(n).  Assumes perfectly phase-coherent
        %   returns and a coherent integrator (e.g. Doppler filter bank).
            arguments
                ~
                nPulses (1,1) double {mustBePositive, mustBeInteger}
            end
            G_dB = 10*log10(nPulses);
        end

        % --------------------------------------------------------------
        function G_dB = noncoherentGain(~, Pfa, Pd, nPulses)
        % noncoherentGain(Pfa, Pd, nPulses)  ->  NCI gain [dB].
        %
        %   Integration gain implied by the Albersheim (1981) approximation:
        %
        %     G_NCI = albersheim(Pfa, Pd, 1) - albersheim(Pfa, Pd, n)
        %
        %   The difference between the n=1 and n-pulse required SNRs is the
        %   effective gain delivered by noncoherent integration.  Integration
        %   loss is coherentGain(n) - noncoherentGain(Pfa, Pd, n).
        %
        %   Valid ranges: see albersheim.
            arguments
                ~
                Pfa     (1,1) double {mustBeInRange(Pfa, 0, 1, "exclusive")}
                Pd      (1,1) double {mustBeInRange(Pd,  0, 1, "exclusive")}
                nPulses (1,1) double {mustBePositive, mustBeInteger}
            end
            G_dB = Radar.RadarReceiver.albersheim(Pfa, Pd, 1) - ...
                   Radar.RadarReceiver.albersheim(Pfa, Pd, nPulses);
        end

        % --------------------------------------------------------------
        function disp(obj)
            fprintf('  RadarReceiver\n');
            fprintf('    NF           = %.4g dB\n',  obj.NF_dB);
            fprintf('    Te           = %.4g K\n',   obj.Te);
            fprintf('    Bn           = %s\n',       Radar.RadarReceiver.fmtFreq(obj.Bn));
            fprintf('    kTBF         = %.4g dBW\n', obj.kTBF_dBW);
            fprintf('    MDS          = %.4g dBW\n', obj.MDS_dBW);
            fprintf('    losses       = %.4g dB\n',  obj.losses_dB);
            fprintf('    T0           = %.4g K\n',   obj.T0);
            if obj.nBits > 0
                fprintf('    nBits        = %d  (DR = %.4g dB)\n', ...
                    obj.nBits, obj.dynamicRange_dB);
            else
                fprintf('    nBits        = ideal\n');
            end
        end

    end   % public methods

    % ==================================================================
    methods (Static)

        % --------------------------------------------------------------
        function SNR1_dB = albersheim(Pfa, Pd, nPulses)
        % albersheim(Pfa, Pd, nPulses)  ->  required per-pulse SNR [dB].
        %
        %   Albersheim (1981) closed-form approximation for noncoherent
        %   square-law integration of nPulses i.i.d. non-fluctuating
        %   (Swerling 0 / Marcum) pulses.
        %
        %   Accuracy: < 0.3 dB for
        %     1e-7 <= Pfa  <= 1e-3
        %     0.1  <= Pd   <= 0.9999
        %     1    <= n    <= 8096
        %
        %   Reference: W. H. Albersheim, "A closed-form approximation to
        %   Robertson's detection characteristics," Proc. IEEE, 1981.
            arguments
                Pfa     (1,1) double {mustBeInRange(Pfa, 0, 1, "exclusive")}
                Pd      (1,1) double {mustBeInRange(Pd,  0, 1, "exclusive")}
                nPulses (1,1) double {mustBePositive, mustBeInteger} = 1
            end
            A = log(0.62 / Pfa);
            B = log(Pd   / (1 - Pd));
            SNR1_dB = -5*log10(nPulses) + ...
                (6.2 + 4.54/sqrt(nPulses + 0.44)) * ...
                log10(A + 0.12*A*B + 1.7*B);
        end

        % --------------------------------------------------------------
        function NF_total_dB = friisChain(NF_dB_vec, gain_dB_vec)
        % friisChain(NF_dB_vec, gain_dB_vec)  ->  cascaded NF [dB].
        %
        %   Computes the Friis noise figure for a K-stage receive chain.
        %
        %   NF_dB_vec   : 1xK  per-stage noise figures [dB]
        %   gain_dB_vec : 1xK  per-stage available gain [dB]
        %                      (passive loss = negative gain)
        %
        %   F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
        %
        %   Reference: H. T. Friis, "Noise figures of radio receivers,"
        %   Proc. IRE, 1944.
            arguments
                NF_dB_vec   (1,:) double
                gain_dB_vec (1,:) double
            end
            assert(numel(NF_dB_vec) == numel(gain_dB_vec), ...
                'RadarReceiver.friisChain: NF_dB_vec and gain_dB_vec must have equal length.');

            F    = 10.^(NF_dB_vec   / 10);
            G    = 10.^(gain_dB_vec / 10);
            Gcum = [1, cumprod(G(1:end-1))];   % cumulative gain before each stage

            NF_total_dB = 10*log10(sum((F - 1) ./ Gcum) + 1);
        end

    end   % static public methods

    % ==================================================================
    methods (Static, Access = private)

        function str = fmtFreq(f)
            if     f >= 1e9, str = sprintf('%.4g GHz', f/1e9);
            elseif f >= 1e6, str = sprintf('%.4g MHz', f/1e6);
            elseif f >= 1e3, str = sprintf('%.4g kHz', f/1e3);
            else,            str = sprintf('%.4g Hz',  f);
            end
        end

    end

end
