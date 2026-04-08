function ula_hpbw(N_max, d_over_lambda, theta0_deg)
    %ULA_HPBW Summary of this function goes here
    %   Detailed explanation goes here
    %
    % INPUT
    %
    %
    % OUTPUT
    %
    %
    % USAGE
    %
    %

    % Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2026-04-07 18:16

    % Exact and approximate HPBW for ULA, broadside or steered.
    % N_max        : compute for N = 2:N_max
    % d_over_lambda: element spacing in wavelengths
    % theta0_deg   : steering angle from broadside axis (0 = broadside)

    if nargin < 1, N_max        = 10;  end
    if nargin < 2, d_over_lambda = 0.5; end
    if nargin < 3, theta0_deg   = 90;  end  % 90 deg = broadside in theta-convention

    kd   = 2*pi * d_over_lambda;
    beta = kd * cosd(theta0_deg);          % progressive phase for steering

    Nvec = (2:N_max)';
    psi_hp_exact = zeros(size(Nvec));
    psi_hp_L1    = zeros(size(Nvec));      % denominator linearized only

    %% ---- sinc root (Level 1 denominator linearization) ----------------------
    % sin(x)/x = 1/sqrt(2), solve once
    x0 = fzero(@(x) sin(x)./x - 1/sqrt(2), 1.2);   % x0 ~ 1.3916

    for k = 1:numel(Nvec)
        N = Nvec(k);

        % --- exact: solve sin(N*u)/sin(u) = N/sqrt(2) for u in (0, pi/N) ----
        % equivalently f(u) = sqrt(2)*sin(N*u) - N*sin(u) = 0
        f = @(u) sqrt(2)*sin(N*u) - N*sin(u);

        % bracket: u=0 is f=0 (trivial), first non-trivial root < pi/N
        % scan for sign change
        u_test = linspace(1e-6, pi/N - 1e-6, 2000);
        fv     = f(u_test);
        idx    = find(diff(sign(fv)) < 0, 1);          % first negative crossing
        u_root = fzero(f, [u_test(idx), u_test(idx+1)]);
        psi_hp_exact(k) = 2*u_root;

        % Level 1
        psi_hp_L1(k) = 2*x0 / N;
    end

    %% ---- map psi_hp -> HPBW in theta-space ----------------------------------
    % psi = kd*cos(theta) - beta  =>  cos(theta) = (psi + beta)/kd
    % theta_hp satisfies psi(theta_hp) = +/- psi_hp around psi=0

    hpbw_exact = hpbw_from_psi(psi_hp_exact, kd, beta);
    hpbw_L1    = hpbw_from_psi(psi_hp_L1,    kd, beta);
    % Level 2: full linearization
    hpbw_L2    = rad2deg(0.886 * 1 ./ (Nvec * d_over_lambda));

    %% ---- table ---------------------------------------------------------------
    fprintf('\nd/lambda = %.3f,  theta0 = %.1f deg\n\n', d_over_lambda, theta0_deg);
    fprintf('%4s  %10s  %10s  %10s  %10s  %8s\n', ...
        'N', 'psi_hp', 'HPBW_exact', 'HPBW_L1', 'HPBW_L2', 'Err_L2%');
    fprintf('%s\n', repmat('-',1,62));
    for k = 1:numel(Nvec)
        err = 100*(hpbw_L2(k) - hpbw_exact(k)) / hpbw_exact(k);
        fprintf('%4d  %10.4f  %10.4f  %10.4f  %10.4f  %8.2f\n', ...
            Nvec(k), psi_hp_exact(k), hpbw_exact(k), hpbw_L1(k), hpbw_L2(k), err);
    end

    %% ---- plot ----------------------------------------------------------------
    figure('Position',[100 100 1100 480]);

    subplot(1,2,1)
    plot(Nvec, hpbw_exact, 'ko-', 'LineWidth',1.5, 'MarkerFaceColor','k', 'DisplayName','Exact'); hold on
    plot(Nvec, hpbw_L1,    'b^--','LineWidth',1.2, 'DisplayName','L1: denom linearized');
    plot(Nvec, hpbw_L2,    'rs:' ,'LineWidth',1.2, 'DisplayName','L2: fully linearized');
    xlabel('N');  ylabel('HPBW (deg)');
    title(sprintf('HPBW vs N  (d/\\lambda=%.2f, \\theta_0=%.0f°)', d_over_lambda, theta0_deg));
    legend('Location','northeast');  grid on;

    subplot(1,2,2)
    err_L1 = 100*(hpbw_L1 - hpbw_exact) ./ hpbw_exact;
    err_L2 = 100*(hpbw_L2 - hpbw_exact) ./ hpbw_exact;
    plot(Nvec, err_L1, 'b^--','LineWidth',1.5,'DisplayName','L1 error'); hold on
    plot(Nvec, err_L2, 'rs:' ,'LineWidth',1.5,'DisplayName','L2 error');
    yline(0,'k-','LineWidth',0.8);
    xlabel('N');  ylabel('Relative error (%)');
    title('Approximation error vs exact');
    legend('Location','southeast');  grid on;

    %% ---- pattern plots for small N ------------------------------------------
    figure('Position',[100 600 1100 380]);
    theta = linspace(0, 180, 3601);
    psi   = kd*cosd(theta) - beta;
    cmap  = lines(min(N_max-1, 6));
    hold on;
    for k = 1:min(numel(Nvec), 6)
        N  = Nvec(k);
        AF = abs(sin(N*psi/2) ./ (N*sin(psi/2) + eps));
        plot(theta, 20*log10(AF + eps), 'Color', cmap(k,:), 'LineWidth',1.2, ...
            'DisplayName', sprintf('N=%d', N));
    end
    yline(-3, 'k--', '-3 dB', 'LineWidth', 0.8);
    xlim([0 180]);  ylim([-40 1]);
    xlabel('\theta (deg)');  ylabel('Normalized AF (dB)');
    title('Array factor magnitude');
    legend('Location','north','NumColumns',3);  grid on;

end

%% =========================================================================
function hpbw = hpbw_from_psi(psi_hp, kd, beta)
    % Exact angle-space HPBW from psi_hp, handling arccos domain clamp.
    % Beam center at psi=0: cos(theta0) = beta/kd
    % Half-power: cos(theta_hp) = (+-psi_hp + beta)/kd

    c_plus  = ( psi_hp + beta) / kd;
    c_minus = (-psi_hp + beta) / kd;

    % clamp to valid arccos domain (scan loss / endfire boundary)
    c_plus  = max(-1, min(1, c_plus));
    c_minus = max(-1, min(1, c_minus));

    theta_plus  = acosd(c_plus);
    theta_minus = acosd(c_minus);

    hpbw = abs(theta_minus - theta_plus);
end


