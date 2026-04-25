classdef RadarReceiver < handle
% RadarReceiver  Facade composing ReceiverNoise, ReceiverChain, and
%               ReceiverDetection into a single radar receiver object.
%
%   Owns the range equation (rangeSNR, maxRange) which requires both the
%   noise model and transmitter parameters.  All other functionality is
%   delegated to the sub-objects exposed as public properties.
%
%   CONSTRUCTION
%     rx = RadarReceiver(noise, chain)
%     rx = RadarReceiver(noise)              – chain defaults applied
%
%   Arguments
%     noise  : ReceiverNoise   object  (required)
%     chain  : ReceiverChain   object  (optional; default constructed if absent)
%
%   SUB-OBJECT ACCESS
%     rx.noise    → ReceiverNoise     (NF, Te, kTBF, MDS, noisePower, snr, friisChain)
%     rx.chain    → ReceiverChain     (process, plot, dynamicRange_dB)
%     rx.detect   → ReceiverDetection (albersheim, coherentGain, noncoherentGain)
%
%   RANGE EQUATION
%     [SNR_dB, Pr_dBW] = rx.rangeSNR(tx, Gt_dBi, Gr_dBi, sigma_m2, R_m)
%     R_m              = rx.maxRange(tx, Gt_dBi, Gr_dBi, sigma_m2, SNRmin_dB)
%
%   PROCESSING CHAIN  (convenience wrappers — pass NF_dB and T0 to chain)
%     C  = rx.process(s_rf, fan)
%     C  = rx.process(s_rf, fan, h_ref=h)
%     rx.plot(C)
%
%   EXAMPLE
%     B   = 10e6;  Tp = 20e-6;  fIF = 60e6;  fan = 600e6;
%
%     noise = ReceiverNoise(4, B, losses_dB=1.5);
%     chain = ReceiverChain(nBits=12, fIF=fIF, fadc=150e6, fbb=25e6);
%     rx    = RadarReceiver(noise, chain);
%     disp(rx)
%
%     % Detection
%     SNR1  = rx.detect.albersheim(1e-6, 0.9, 1)
%     G_nci = rx.detect.noncoherentGain(1e-6, 0.9, 16)
%
%     % Processing chain
%     t    = (0:round(Tp*fan)-1).' / fan;
%     s_rf = cos(2*pi*fIF*t + pi*(B/Tp)*(t - Tp/2).^2);
%     t_r  = (0:round(Tp*chain.fbb)-1).' / chain.fbb;
%     h    = exp(1j*pi*(B/Tp)*(t_r - Tp/2).^2);
%     C    = rx.process(s_rf, fan, h_ref=h);
%     rx.plot(C);

    % ------------------------------------------------------------------
    properties (SetAccess = private)
        noise    % ReceiverNoise
        chain    % ReceiverChain
        detect   % ReceiverDetection  (stateless; shared instance)
    end

    properties (Constant, Access = private)
        c0 = 299792458;
    end

    % ==================================================================
    methods

        function obj = RadarReceiver(noise, chain)
            arguments
                noise           Radar.ReceiverNoise
                chain           Radar.ReceiverChain  = Radar.ReceiverChain()
            end
            obj.noise  = noise;
            obj.chain  = chain;
            obj.detect = Radar.ReceiverDetection;
        end

        % ── Range equation ────────────────────────────────────────────

        function [SNR_dB, Pr_dBW] = rangeSNR(obj, tx, Gt_dBi, Gr_dBi, sigma_m2, R_m)
        % rangeSNR(tx, Gt_dBi, Gr_dBi, sigma_m2, R_m)
        %
        %   Monostatic radar range equation at range(s) R_m [m].
        %   System loss = tx.losses_dB + noise.losses_dB.
        %   Includes waveform compression gain (TBP) from transmitter.
        %
        %   Returns SNR_dB and received power Pr_dBW, both Nx1.
            arguments
                obj
                tx              Radar.RadarTransmitter
                Gt_dBi    (1,1) double
                Gr_dBi    (1,1) double
                sigma_m2  (1,1) double {mustBePositive}
                R_m       (:,1) double {mustBePositive}
            end
            Gt   = 10^(Gt_dBi / 10);
            Gr   = 10^(Gr_dBi / 10);
            Ltot = 10^((tx.losses_dB + obj.noise.losses_dB) / 10);
            lam  = tx.lambda;
            TBP  = tx.timeBandwidthProduct;

            Pr_W   = (tx.Pt * Gt * Gr * lam^2 * sigma_m2 * TBP) ./ ...
                     ((4*pi)^3 * R_m.^4 * Ltot);
            Pr_dBW = 10*log10(Pr_W);
            SNR_dB = Pr_dBW - obj.noise.kTBF_dBW;
        end

        function R_m = maxRange(obj, tx, Gt_dBi, Gr_dBi, sigma_m2, SNRmin_dB)
        % maxRange(tx, Gt_dBi, Gr_dBi, sigma_m2, SNRmin_dB)  ->  [m]
        %
        %   Closed-form detection range for a required single-pulse SNR.
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
            Ltot   = 10^((tx.losses_dB + obj.noise.losses_dB) / 10);

            num = tx.Pt * Gt * Gr * tx.lambda^2 * sigma_m2;
            den = (4*pi)^3 * obj.noise.kTBF_W * Ltot * SNRmin;
            R_m = (num / den)^0.25;
        end

        % ── Chain convenience wrappers ────────────────────────────────

        function C = process(obj, s_rf, fan, opts)
        % process(s_rf, fan)           — delegates to chain.process,
        % process(s_rf, fan, h_ref=h)    injecting NF_dB and T0 from noise.
            arguments
                obj
                s_rf  (:,1) double
                fan   (1,1) double {mustBePositive}
                opts.h_ref (:,1) {mustBeNumeric} = []
            end
            C = obj.chain.process(s_rf, fan, ...
                NF_dB  = obj.noise.NF_dB, ...
                T0     = obj.noise.T0,    ...
                h_ref  = opts.h_ref);
        end

        function plot(obj, C)
        % plot(C)  — delegates to chain.plot.
            obj.chain.plot(C);
        end

        % ── Display ───────────────────────────────────────────────────

        function disp(obj)
            fprintf('  RadarReceiver\n');
            fprintf('  ── Noise chain ─────────────────────────────\n');
            disp(obj.noise);
            fprintf('  ── Processing chain ────────────────────────\n');
            disp(obj.chain);
        end

    end

end
