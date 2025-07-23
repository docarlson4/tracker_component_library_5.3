clc
clear
close all

cd(fileparts(which(mfilename)))

startup_tcl

dbstop if error

% Doug Carlson (doug.o.carlson@gmail.com), 2024-04-19 16:30

%% Parameters
rng(0)

% Number of scans
num_scans = 10000;

% Averaging length
L = 1000;

% Region of interest
Lx = 12000;

% Clutter density
lambda_lo = 1e-3; % Nc/m^2
lambda_hi = 2e-3; % Nc/m^2

mmt_space_hi = [ -Lx/2; -Lx/4];
mmt_space_lo = [ -Lx/4;  Lx/2];

lx = 1000;

%% Clutter Map Region

xc = 0:lx:Lx;
xc = xc(2:end) - diff(xc)/2;
xc = xc - mean(xc);

num_cells = numel(xc);

mmt_vol_lo = prod(diff(mmt_space_lo)); % m^2
mmt_vol_hi = prod(diff(mmt_space_hi)); % m^2

zClut = cell(num_scans, 1);
for kS = 1:num_scans
    % Determine the number of false alarms to generate.
    num_lo = PoissonD.rand(1, mmt_vol_lo * lambda_lo);
    num_hi = PoissonD.rand(1, mmt_vol_hi * lambda_hi);

    x_lo = UniformD.rand(num_lo, mmt_space_lo(:,1));
    x_hi = UniformD.rand(num_hi, mmt_space_hi(:,1));

    zClut{kS} = [x_lo, x_hi];

end

%% Classic Clutter Map
map_classic = zeros(1, num_cells);
for kS = 1:num_scans
    idx = min(kS, L);
    alp = (idx-1)/idx; % Compute alpha outside the loop
    lgcl_x = abs(bsxfun(@minus, zClut{kS}(1,:), xc(:))) < lx/2;
    mu = sum(lgcl_x, 2)'; % Sum over the second dimension
    mu(mu == 0) = 1/lx;
    map_classic = alp .* map_classic + (1-alp) .* mu;
end
map_classic = map_classic/(lx);

%% Spatial Clutter Map

map_spatial = zeros(1, num_cells);
mu = zeros(L,num_cells);
for kS = 1:num_scans

    xm = zClut{kS}(1,:);

    for kC = 1:num_cells
        [del_x, idx_x] = min( abs(xm - xc(kC)) );

        rhs = xc(kC) + del_x;
        if xc(kC) == max(xc)
            rhs = min(xc(kC) + del_x, xc(kC) + lx/2);
        end
        
        lhs = xc(kC) - del_x;
        if xc(kC) == min(xc)
            lhs = max(xc(kC) - del_x, xc(kC) - lx/2);
        end

        if isempty(lhs) && isempty(rhs)
            mu(min(kS,L), kC) = 0;
        else
            mu(min(kS,L), kC) = rhs - lhs;
        end
    end
    if kS >= L
        map_spatial = mean(mu, 1);
        mu = circshift(mu, [-1,0]);
    end
end

%% Plots

figure('Position',[200,200,560,420])
hold on, grid on, box on
scatter(xc,map_classic, DisplayName="Classic")
scatter(xc,1./map_spatial,'*', DisplayName="Spatial")
legend
ylim([0,sum([lambda_hi, lambda_lo])])