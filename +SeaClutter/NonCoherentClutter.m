classdef NonCoherentClutter < handle
    %NONCOHERENTCLUTTER This first level corresponds to the generation of
    %sea clutter data intended for non-Doppler processing using a single
    %channel antenna: thus, the CPI is equal to one single PRI (NFFT = 1).
    %
    % CONSTRUCTION
    %   a = NonCoherentClutter()
    %   a = NonCoherentClutter(Name=Value, ...)
    % 
    % NAMED PARAMETERS (with defaults)
    %   TextureModel (1,1) string   - ["Gamma", "InverseGamma"]
    %   ShapeParam (1,1) double     - [nu]
    %   CorrCoeffRng (1,1) function_handle = @deal
    % 
    % 
    % 
    % 
    % 
    % [0] Kemkemian, S., Lupinski, L., Corretja, V., Cottron, R., & Watts,
    % S. (2015, May). Performance assessment of multi-channel radars using
    % simulated sea clutter. In 2015 IEEE Radar Conference (RadarCon) (pp.
    % 1015-1020). IEEE.

    % -----------------------------------------------------------------
    %  Public properties
    % -----------------------------------------------------------------
    properties
        TextureModel %(1,1) string = "Gamma" % ["Gamma", "InverseGamma"]
        ShapeParam %(1,1) double = 1    % [nu]
        CorrCoeffRng %(1,1) function_handle = @deal
        NumRange %(1,1) double {mustBePositive} = 64;
        NumAzimuth %(1,1) double {mustBePositive} = 128;
    end

    % -----------------------------------------------------------------
    %  Private properties
    % -----------------------------------------------------------------
    properties (Access = private)
    end

    % -----------------------------------------------------------------
    %  Dependent property — recomputed from Wavelength automatically
    % -----------------------------------------------------------------
    properties (Dependent, GetAccess = public, SetAccess = private)
        ScaleParam
        Data
    end

    methods
        function obj = NonCoherentClutter(opts)
            %NONCOHERENTCLUTTER Construct an instance of this class
            arguments
                opts.TextureModel (1,1) string = "Gamma" % ["Gamma", "InverseGamma"]
                opts.ShapeParam (1,1) double {mustBePositive} = 1.0 % [nu]
                opts.CorrCoeffRng (1,1) function_handle = @deal
                opts.NumRange (1,1) double {mustBePositive} = 64;
                opts.NumAzimuth (1,1) double {mustBePositive} = 128;
            end
            fields = fieldnames(opts);
            for i = 1:numel(fields)
                obj.(fields{i}) = opts.(fields{i});
            end

        end

        % ---- Dependent getter --------------------------------------- %
        function value = get.ScaleParam(obj)
            value = 1/obj.ShapeParam;
        end

        function value = get.Data(obj)
            Nr = obj.NumRange;
            Naz = obj.NumAzimuth;
            x = complex(randn(Nr,Naz), randn(Nr,Naz))/sqrt(2);

            rho_x = obj.CorrCoeffRng(0:Nr-1)' * ones(1,Naz);

            params.type_out = obj.TextureModel;
            params.alpha = obj.ShapeParam;  % shape parameter
            params.beta = obj.ScaleParam; % scale parameter
            num_gauss = 300;
            MC_INTEG = 0; % perform Monte Carlo integration or Hermite
            [lut_rho, lut_RG] = SeaClutter.CorrCoefMapping(params, num_gauss, MC_INTEG);

            % Pre-warp desired ACF (clamp to achievable range [0, max(lut_rho)])
            rho_x_cl = min(max(rho_x, 0), lut_rho(end));
            RG_x = reshape(interp1(lut_rho, lut_RG, rho_x_cl(:), ...
                'pchip', 0), Nr, Naz);
            RG_x(1,1) = 1;                  % zero-lag must be 1

            % Generate 2D correlated Gaussian field with ACF = RG_x
            g_x = gen_corr_gaussian_2d(RG_x, Nr, Naz);

            % MNLT: N(0,1) -> Gamma(nu, 1/nu)  [unit mean, variance = 1/nu]
            u_x    = min(max(GaussianD.CDF(g_x), 1e-15), 1-1e-15);
            tau = GammaD.invCDF(u_x, obj.ShapeParam, obj.ScaleParam);   % unit-mean texture

            value = sqrt(tau) .* x;
        end

        function PlotClutterIntensity(obj)
            Nr = obj.NumRange;
            Naz = obj.NumAzimuth;
            figure('Name','Range-Time Intensity','Position',[50 150 900 400]);
            imagesc((1:Naz), (1:Nr), 10*log10(abs(obj.Data).^2 + 1e-12));
            colorbar; colormap(jet);
            xlabel('Azimuth Bin'); ylabel('Range gate');
            title(sprintf('Clutter intensity [dB]  \\nu=%.1f, CNR=%.0f dB', ...
                obj.ShapeParam, nan));
            set(gca,'YDir','normal');
        end

        function VerifyTextureACF(obj)
            Nr = obj.NumRange;
            Naz = obj.NumAzimuth;
            figure('Name','Texture ACF','Position',[550 150 480 350]);
            x_unit    = obj.Data;
            acf_emp_r = zeros(2*Nr-1,1);
            for kAz = 1:Naz
                acf_emp_r = acf_emp_r + xcorr(x_unit(:,kAz));
            end
            lags_r = 0:Nr/2-1;
            acf_emp_r = acf_emp_r/Nr;
            plot(lags_r, acf_emp_r(1:Nr/2), 'b.-', 'MarkerSize', 8); hold on;
            plot(lags_r, obj.CorrCoeffRng(lags_r), 'r--', 'LineWidth', 2);
            xlabel('Range lag (cells)'); ylabel('Normalised ACF \rho_y');
            legend('Simulated', 'Target Gaussian ACF', 'Location', 'NE');
            title('Texture range ACF');
        end

    end % public methods

    methods (Access = private)
    end % private methods

end % classdef

%% =========================================================================
%  LOCAL FUNCTION: gen_corr_gaussian_2d
%% =========================================================================
function G = gen_corr_gaussian_2d(R2d_circ, Nr, Np)
%GEN_CORR_GAUSSIAN_2D  2-D zero-mean unit-variance real Gaussian field.
%
% R2d_circ: [Nr x Np] circular ACF (element (1,1) is zero-lag = 1),
%           arranged for MATLAB's fft2 (wrap-around lags).
%
% Method: spectral shaping of white Gaussian noise.
%   G = real(ifft2(sqrt(S) .* fft2(Z)))
% where S = fft2(R2d_circ) is the 2-D PSD.
%
% Variance derivation: Var(G) = (1/(NrNp)^2) * sum_{k,l} S(k,l)*NrNp = 1,
% confirming unit variance exactly for the circular-stationary model.

S = max(real(fft2(R2d_circ)), 0);   % 2-D PSD; zero negative numerical noise
Z = randn(Nr, Np);
G = real(ifft2(sqrt(S) .* fft2(Z)));

% Enforce empirical N(0,1) to remove finite-sample bias (optional but robust)
G = (G - mean(G(:))) / std(G(:));
end

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2026-03-20 13:40
