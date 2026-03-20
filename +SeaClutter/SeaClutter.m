classdef SeaClutter < handle
    properties(Access = public)
        TexturePDF (1,1) string = "Gamma"
        Data (:,:) double = []
        Shape (1,1) double {mustBePositive} = 1.0
        NumRange (1,1) double {mustBePositive} = 64;
        NumAzimuth (1,1) double {mustBePositive} = 128;
        Enabled (1,1) logical = true
    end

    properties (Dependent, GetAccess = public, SetAccess = private)
        Scale
    end

    properties (Dependent)
        DefaultDataDir
    end

    methods
        function p = get.DefaultDataDir(~)
            if ispc
                p = "c:\Users\dougo\source\repos\tracking\tracker_component_library_5.3\";
            else
                p = "/usr/local/share/data";
            end
        end
    end

    methods
        function obj = SeaClutter(opts)
            % Define valid arguments dynamically from class properties
            arguments
                opts.TexturePDF (1,1) string = "Gamma"
                opts.Data (:,:) double = []
                opts.Shape (1,1) double {mustBePositive} = 1.0
                opts.NumRange (1,1) double {mustBePositive} = 64;
                opts.NumAzimuth (1,1) double {mustBePositive} = 128;
                opts.Enabled (1,1) logical = true
            end

            % Assign all provided name-value pairs
            propNames = fieldnames(opts);
            for i = 1:numel(propNames)
                obj.(propNames{i}) = opts.(propNames{i});
            end

            % Add to MATLAB path if not already there
            if ~contains(path, obj.DefaultDataDir)
                addpath(genpath(obj.DefaultDataDir));
            end
        end
        function value = get.Scale(obj)
            value = 1/obj.Shape;
        end

        function PlotTexture(obj, type)
            arguments
                obj
                type string = "Histogram"
            end
            figure
            hold on, grid on, box on
            switch type
                case "Histogram"
                    nu = obj.Shape;
                    b = obj.Scale;
                    yr = GammaD.rand([obj.NumRange,obj.NumAzimuth], nu, b);
                    histogram(yr, Normalization="pdf")
                    y = 0:0.01:max(yr(:));
                    plot(y, obj.gammaPDF(y))
                case "TimeSeries"
                case "Image"
                    imagesc(dB(obj.Data.^2, -100))
                    axis image
                    colorbar
                    colormap("abyss")
            end
        end
        function data = gammaPDF(obj, y)
            nu = obj.Shape;
            b = obj.Scale;
            data = b.^-nu./gamma(nu) .* y.^(nu-1) .* exp(-y/b);
        end

        function obj = GenClutter(obj)
            x = synth2D_ACF(obj);
            naz = obj.NumAzimuth;
            nrs = obj.NumRange;
            tau = GammaD.rand([nrs,naz], obj.Shape, obj.Scale);
            obj.Data = sqrt(tau) .* x;
        end
    end

    methods(Access=private)
        function f = synth2D_ACF_old(obj)
            naz = obj.NumAzimuth;
            nrs = obj.NumRange;

            m = -nrs/2:nrs/2-1;
            n = -naz/2:naz/2-1;
            rho_rs = exp(-abs(m)/8).*cos(pi*m/8);
            rho_az = exp(-abs(n)/40);
            R = rho_rs' * rho_az;
            [Nx,Ny] = size(R);
            S = real(fft2(R)); % PSD

            % Ensure non-negative PSD
            S = max(S,0);

            % Random phases with Hermitian symmetry
            phase = randn(Nx,Ny) + 1i*randn(Nx,Ny);
            phase(1,1) = real(phase(1,1)); % DC term real
            for ix = 2:Nx
                for iy = 2:Ny
                    kx_sym = mod(Nx-ix+2,Nx)+1;
                    ky_sym = mod(Ny-iy+2,Ny)+1;
                    phase(kx_sym,ky_sym) = conj(phase(ix,iy));
                end
            end

            F = sqrt(S).*phase;
            f = real(ifft2(F));
        end

        function [field, acf_est] = synth2D_ACF(obj)
            %(rho2D, dx, dy)
            % synth2D_ACF - Synthesize a real 2D random field with a target 2D ACF
            %
            % Inputs:
            %   rho2D : Nx-by-Ny matrix containing the target 2D autocorrelation function
            %   dx    : spatial resolution in x
            %   dy    : spatial resolution in y
            %
            % Outputs:
            %   field   : Nx-by-Ny real-valued synthesized field
            %   acf_est : Nx-by-Ny empirical autocorrelation of the field
            %
            % Example:
            %   Nx = 128; Ny = 128;
            %   dx = 0.1; dy = 0.1;
            %   [X,Y] = meshgrid((0:Nx-1)*dx,(0:Ny-1)*dy);
            %   rho2D = exp(-sqrt((X-mean(X(:))).^2 + (Y-mean(Y(:))).^2)/2);
            %   [field,acf_est] = synth2D_ACF(rho2D, dx, dy);

            naz = obj.NumAzimuth;
            nrs = obj.NumRange;

            m = -nrs/2:nrs/2-1;
            n = -naz/2:naz/2-1;
            rho_rs = exp(-abs(m)/8).*cos(pi*m/8);
            rho_az = exp(-abs(n)/20);
            rho2D = rho_rs' * rho_az;
            
            [Nx, Ny] = size(rho2D);
            dx = 1;
            dy = 1;

            % Step 1 — Compute 2D PSD from target ACF
            S_two_sided = real(fft2(ifftshift(rho2D))) * dx * dy;
            S_two_sided = max(S_two_sided,0); % numerical safety

            % Step 2 — Generate Hermitian random phases
            rand_phase = randn(Nx,Ny) + 1i*randn(Nx,Ny);

            for ix = 1:Nx
                for iy = 1:Ny
                    ix_sym = mod(Nx-ix+2,Nx)+1;
                    iy_sym = mod(Ny-iy+2,Ny)+1;
                    rand_phase(ix_sym,iy_sym) = conj(rand_phase(ix,iy));
                end
            end

            % Step 3 — Apply magnitude sqrt(PSD)
            F = sqrt(S_two_sided) .* rand_phase;

            % Step 4 — Force DC and Nyquist bins to be real
            F(1,1) = real(F(1,1));
            if mod(Nx,2)==0
                F(Nx/2+1,:) = real(F(Nx/2+1,:));
            end
            if mod(Ny,2)==0
                F(:,Ny/2+1) = real(F(:,Ny/2+1));
            end

            % Step 5 — Inverse FFT to get real-valued field
            field = real(ifft2(F));

            % Step 6 — Empirical autocorrelation
            acf_est = fftshift(ifft2(abs(fft2(field)).^2) * dx * dy);

        end

        function x = genSIRP(obj)
            naz = obj.NumAzimuth;
            nhaz = naz/2+1;
            nrs = obj.NumRange;
            nhrs = nrs/2+1;

            m = -nrs/2:nrs/2-1;
            n = -naz/2:naz/2-1;
            rho_rs = exp(-abs(m)/8);
            rho_az = exp(-abs(n)/4);

            S_rs_2s = real(fft(ifftshift(rho_rs)));
            S_az_2s = real(fft(ifftshift(rho_az)));

            S_rs_pos = S_rs_2s(1:nhrs);
            S_az_pos = S_az_2s(1:nhaz);

            X_rs_pos = zeros(1,nhrs);
            X_rs_pos(1) = sqrt(S_rs_pos(1)/nrs) * randn();
            X_rs_pos(end) = sqrt(S_rs_pos(end)/nrs) * randn();
            k = 2:(nhrs-1);
            sigma = sqrt( S_rs_pos(k) / nrs / 2 );   % /2 because complex has two quadratures
            X_rs_pos(k) = sigma .* ( randn(size(k)) + 1i*randn(size(k)) );
            X_rs_full = [ X_rs_pos, conj( fliplr( X_rs_pos(2:nhrs-1) ) ) ];

            X_az_pos = zeros(1,nhaz);
            X_az_pos(1) = sqrt(S_az_pos(1)/naz) * randn();
            X_az_pos(end) = sqrt(S_az_pos(end)/naz) * randn();
            k = 2:(nhaz-1);
            sigma = sqrt( S_az_pos(k) / naz / 2 );   % /2 because complex has two quadratures
            X_az_pos(k) = sigma .* ( randn(size(k)) + 1i*randn(size(k)) );
            X_az_full = [ X_az_pos, conj( fliplr( X_az_pos(2:nhaz-1) ) ) ];

            % x = complex(randn(naz,naz), randn(naz,naz))/sqrt(2);
            X_full = X_rs_full.' * X_az_full;
            x = real(ifft2( X_full ));      % should be real by construction

        end
    end
end
