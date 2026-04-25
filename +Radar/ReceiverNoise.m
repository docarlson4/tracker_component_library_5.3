classdef ReceiverNoise < handle
% ReceiverNoise  Receiver noise chain: NF, thermal noise, sensitivity.
%
%   CONSTRUCTION
%     n = ReceiverNoise(NF_dB, Bn)
%     n = ReceiverNoise(NF_dB, Bn, losses_dB=L, T0=T)
%
%   Required
%     NF_dB   – receiver noise figure  [dB]
%     Bn      – noise bandwidth        [Hz]
%
%   Optional
%     losses_dB – ohmic Rx losses      [dB]   (default: 0)
%     T0        – reference temp       [K]    (default: 290)
%
%   STATIC
%     ReceiverNoise.friisChain(NF_dB_vec, gain_dB_vec)  ->  cascaded NF [dB]
%
%   EXAMPLE
%     n = ReceiverNoise(4, 5e6, losses_dB=1.5);
%     disp(n)
%     NF_cas = ReceiverNoise.friisChain([2 8 5], [20 -6 30])

    properties (SetAccess = private)
        NF_dB       % Noise figure          [dB]
        Bn          % Noise bandwidth       [Hz]
        losses_dB   % Ohmic Rx losses       [dB]
        T0          % Reference temperature [K]
    end

    properties (Dependent)
        F           % Noise figure (linear)
        Te          % Effective noise temperature  T0*(F-1)  [K]
        kTBF_W      % Thermal noise power at Bn  [W]
        kTBF_dBW    % Thermal noise power at Bn  [dBW]
        MDS_dBW     % Min detectable signal (SNR=0 dB + losses)  [dBW]
        rangeRes    % Range resolution  c0/(2*Bn)  [m]
    end

    properties (Constant, Access = private)
        kB = 1.380649e-23;
        c0 = 299792458;
    end

    methods

        function obj = ReceiverNoise(NF_dB, Bn, opts)
            arguments
                NF_dB             (1,1) double {mustBeNonnegative}
                Bn                (1,1) double {mustBePositive}
                opts.losses_dB    (1,1) double {mustBeNonnegative} = 0
                opts.T0           (1,1) double {mustBePositive}    = 290
            end
            obj.NF_dB     = NF_dB;
            obj.Bn        = Bn;
            obj.losses_dB = opts.losses_dB;
            obj.T0        = opts.T0;
        end

        % ── Dependent getters ─────────────────────────────────────────
        function v = get.F(obj),        v = 10^(obj.NF_dB/10);               end
        function v = get.Te(obj),       v = obj.T0*(obj.F - 1);              end
        function v = get.kTBF_W(obj),   v = obj.kB*obj.T0*obj.Bn*obj.F;     end
        function v = get.kTBF_dBW(obj), v = 10*log10(obj.kTBF_W);           end
        function v = get.MDS_dBW(obj),  v = obj.kTBF_dBW + obj.losses_dB;   end
        function v = get.rangeRes(obj), v = obj.c0 / (2*obj.Bn);            end

        % ── Instance methods ──────────────────────────────────────────

        function [Pn_W, Pn_dBW] = noisePower(obj, BW)
        % noisePower(BW)  ->  noise power [W] and [dBW] at bandwidth BW.
            arguments, obj; BW (1,1) double {mustBePositive}; end
            Pn_W   = obj.kB * obj.T0 * BW * obj.F;
            Pn_dBW = 10*log10(Pn_W);
        end

        function [SNR, SNR_dB] = snr(obj, Pr_W)
        % snr(Pr_W)  ->  SNR [linear] and [dB] at noise bandwidth Bn.
        %   Pr_W may be a vector.  Includes Rx losses.
            arguments, obj; Pr_W (:,1) double {mustBePositive}; end
            Lrx    = 10^(obj.losses_dB / 10);
            SNR    = Pr_W ./ (obj.kTBF_W * Lrx);
            SNR_dB = 10*log10(SNR);
        end

        function disp(obj)
            fprintf('  ReceiverNoise\n');
            fprintf('    NF        = %.4g dB\n',  obj.NF_dB);
            fprintf('    Te        = %.4g K\n',   obj.Te);
            fprintf('    Bn        = %s\n',       Radar.ReceiverNoise.fmtHz(obj.Bn));
            fprintf('    kTBF      = %.4g dBW\n', obj.kTBF_dBW);
            fprintf('    MDS       = %.4g dBW\n', obj.MDS_dBW);
            fprintf('    losses    = %.4g dB\n',  obj.losses_dB);
            fprintf('    T0        = %.4g K\n',   obj.T0);
            fprintf('    rangeRes  = %.4g m\n',   obj.rangeRes);
        end

    end

    methods (Static)

        function NF_dB = friisChain(NF_dB_vec, gain_dB_vec)
        % friisChain(NF_dB_vec, gain_dB_vec)  ->  cascaded NF [dB].
        %
        %   F_total = F1 + (F2-1)/G1 + (F3-1)/(G1*G2) + ...
        %
        %   gain_dB_vec: passive loss = negative gain.
            arguments
                NF_dB_vec   (1,:) double
                gain_dB_vec (1,:) double
            end
            assert(numel(NF_dB_vec) == numel(gain_dB_vec), ...
                'ReceiverNoise.friisChain: vectors must be the same length.');
            F    = 10.^(NF_dB_vec   / 10);
            G    = 10.^(gain_dB_vec / 10);
            Gcum = [1, cumprod(G(1:end-1))];
            NF_dB = 10*log10(sum((F-1)./Gcum) + 1);
        end

    end

    methods (Static, Access = private)
        function s = fmtHz(f)
            if     f >= 1e9, s = sprintf('%.4g GHz', f/1e9);
            elseif f >= 1e6, s = sprintf('%.4g MHz', f/1e6);
            elseif f >= 1e3, s = sprintf('%.4g kHz', f/1e3);
            else,            s = sprintf('%.4g Hz',  f);
            end
        end
    end

end
