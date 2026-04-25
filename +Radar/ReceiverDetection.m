classdef ReceiverDetection
% ReceiverDetection  Detection theory utilities: Albersheim approximation,
%                    coherent and noncoherent integration gains.
%
%   All methods are static — no instantiation required.
%
%   METHODS
%     ReceiverDetection.albersheim(Pfa, Pd, n)   ->  required SNR [dB]
%     ReceiverDetection.coherentGain(n)           ->  CI gain [dB]
%     ReceiverDetection.noncoherentGain(Pfa,Pd,n) ->  NCI gain [dB]
%     ReceiverDetection.integrationLoss(Pfa,Pd,n) ->  loss vs CI [dB]
%
%   EXAMPLE
%     Pfa = 1e-6;  Pd = 0.9;  n = 16;
%     snr1   = ReceiverDetection.albersheim(Pfa, Pd, 1)
%     snrN   = ReceiverDetection.albersheim(Pfa, Pd, n)
%     G_ci   = ReceiverDetection.coherentGain(n)
%     G_nci  = ReceiverDetection.noncoherentGain(Pfa, Pd, n)
%     L_int  = ReceiverDetection.integrationLoss(Pfa, Pd, n)

    methods (Static)

        function SNR_dB = albersheim(Pfa, Pd, nPulses)
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
            arguments
                Pfa     (1,1) double {mustBeInRange(Pfa, 0, 1, "exclusive")}
                Pd      (1,1) double {mustBeInRange(Pd,  0, 1, "exclusive")}
                nPulses (1,1) double {mustBePositive, mustBeInteger} = 1
            end
            A = log(0.62 / Pfa);
            B = log(Pd   / (1 - Pd));
            SNR_dB = -5*log10(nPulses) + ...
                (6.2 + 4.54/sqrt(nPulses + 0.44)) * ...
                log10(A + 0.12*A*B + 1.7*B);
        end

        function G_dB = coherentGain(nPulses)
        % coherentGain(nPulses)  ->  CI gain [dB].
        %
        %   G_CI = 10*log10(n).  Assumes phase-coherent returns and a
        %   coherent integrator (e.g. Doppler filter bank).
            arguments
                nPulses (1,1) double {mustBePositive, mustBeInteger}
            end
            G_dB = 10*log10(nPulses);
        end

        function G_dB = noncoherentGain(Pfa, Pd, nPulses)
        % noncoherentGain(Pfa, Pd, nPulses)  ->  NCI gain [dB].
        %
        %   Difference in required SNR between n=1 and n pulses via Albersheim.
            arguments
                Pfa     (1,1) double {mustBeInRange(Pfa, 0, 1, "exclusive")}
                Pd      (1,1) double {mustBeInRange(Pd,  0, 1, "exclusive")}
                nPulses (1,1) double {mustBePositive, mustBeInteger}
            end
            G_dB = Radar.ReceiverDetection.albersheim(Pfa, Pd, 1) - ...
                   Radar.ReceiverDetection.albersheim(Pfa, Pd, nPulses);
        end

        function L_dB = integrationLoss(Pfa, Pd, nPulses)
        % integrationLoss(Pfa, Pd, nPulses)  ->  NCI loss relative to CI [dB].
        %
        %   L = coherentGain(n) - noncoherentGain(Pfa, Pd, n)
        %   Always >= 0 dB.
            arguments
                Pfa     (1,1) double {mustBeInRange(Pfa, 0, 1, "exclusive")}
                Pd      (1,1) double {mustBeInRange(Pd,  0, 1, "exclusive")}
                nPulses (1,1) double {mustBePositive, mustBeInteger}
            end
            L_dB = Radar.ReceiverDetection.coherentGain(nPulses) - ...
                   Radar.ReceiverDetection.noncoherentGain(Pfa, Pd, nPulses);
        end

    end

end
