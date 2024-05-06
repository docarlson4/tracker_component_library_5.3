clc
clear
close all

cd(fileparts(which(mfilename)))

startup_tcl

dbstop if error

% Constants

% Doug Carlson (doug.o.carlson@gmail.com), 2024-04-19 16:30

%% Parameters
rng(0)

% Number of scans
num_scans = 100;

% Averaging length
L = 80;

% Region of interest
Lx = 1200;
Ly = 1200;

% Clutter density
lambda_lo = 1e-5; % Nc/m^2
lambda_hi = 5e-5; % Nc/m^2

mmt_space_hi = [ -Lx/2  Ly/6; Lx/2 Ly/2];
mmt_space_lo = [ -Lx/2 -Ly/2; Lx/2 Ly/6];

lx = 100;
ly = 100;

%% Clutter Map Region

xc = 0:lx:Lx;
xc = xc(2:end) - diff(xc)/2;
xc = xc - mean(xc);
num_x = numel(xc);

yc = 0:ly:Ly;
yc = yc(2:end) - diff(yc)/2;
yc = yc - mean(yc);
num_y = numel(yc);

[Xc,Yc] = meshgrid(xc,yc);
num_cells = numel(Xc);

mmt_vol_lo = prod(diff(mmt_space_lo)); % m^2
mmt_vol_hi = prod(diff(mmt_space_hi)); % m^2

zClut = cell(num_scans, 1);
for kS = 1:num_scans
    % Determine the number of false alarms to generate.
    num_lo = PoissonD.rand(1, mmt_vol_lo * lambda_lo);
    num_hi = PoissonD.rand(1, mmt_vol_hi * lambda_hi);

    x_lo = UniformD.rand(num_lo, mmt_space_lo(:,1));
    y_lo = UniformD.rand(num_lo, mmt_space_lo(:,2));
    x_hi = UniformD.rand(num_hi, mmt_space_hi(:,1));
    y_hi = UniformD.rand(num_hi, mmt_space_hi(:,2));

   zClut{kS} = [x_lo, x_hi; y_lo, y_hi];

end

%% Confirm ClutterEstimation Class
cm_obj = Tracker.ClutterEstimation("Type","Classical", ...
    "AveragingLength", L, "MmtRegion", [-Lx/2, -Ly/2; Lx/2, Ly/2], ...
    "CellSize", [lx,ly])

%% Classic Clutter Map
map_classic = zeros(num_y, num_x);
for kS = 1:num_scans
    idx = min(kS, L);
    alp = (idx-1)/idx; % Compute alpha outside the loop
    lgcl_x = abs(bsxfun(@minus, zClut{kS}(1,:), Xc(:))) < lx/2;
    lgcl_y = abs(bsxfun(@minus, zClut{kS}(2,:), Yc(:))) < ly/2;
    mu = sum(lgcl_y & lgcl_x, 2); % Sum over the second dimension
    mu = reshape(mu, [num_y, num_x]);
    map_classic = alp .* map_classic + (1-alp) .* mu;

    cm_obj.ClutterMap(zClut{kS}, kS);
end
map_classic = map_classic/(lx*ly);

%% Spatial Clutter Map

% rat = lx/ly;
% map_spatial = zeros(num_y, num_x);
% for kS = 1:num_scans
%     idx = min(kS, L);
%     alp = (idx-1)/idx; % Compute alpha outside the loop
% 
%     xm = zClut{kS}(1,:);
%     ym = zClut{kS}(2,:);
% 
%     mu = zeros(num_cells,1);
%     for kC = 1:num_cells
%         [~, idx_x] = min( abs(xm - Xc(kC)) );
%         [~, idx_y] = min( abs(ym - Yc(kC)) );
%         if abs(xm(idx_x) - Xc(kC)) < abs(ym(idx_y) - Yc(kC)) * rat
%             hgt_m = abs(ym(idx_y) - Yc(kC)) * 2;
%             wid_m = hgt_m * rat;
%         else
%             wid_m = abs(xm(idx_x) - Xc(kC)) * 2;
%             hgt_m = wid_m / rat;
%         end
%         mu(kC) = wid_m * hgt_m;
%     end
% 
%     mu = reshape(mu, [num_y, num_x]);
%     map_spatial = alp .* map_spatial + (1-alp) .* mu;
% end

%% Temporal Clutter Map

% map_temporal = zeros(num_y, num_x);
% tau_last_mmt = zeros(num_y * num_x, 1);
% for kS = 1:num_scans
%     idx = min(kS, L);
%     alp = (idx-1)/idx; % Compute alpha outside the loop
%     lgcl_x = abs(bsxfun(@minus, zClut{kS}(1,:), X(:))) < lx/2;
%     lgcl_y = abs(bsxfun(@minus, zClut{kS}(2,:), Y(:))) < ly/2;
%     nu = sum(lgcl_y & lgcl_x, 2); % Sum over the second dimension
%     lgcl_mmt_in_cell = (nu > 0);
%     tau = nu();
% 
%     nu = reshape(nu, [num_y, num_x])
%     tau = reshape(tau, [num_y, num_x])
% 
%     map_classic = alp .* map_classic + (1-alp) .* mu;
% end

%% Plots

figure('Position',[100,200,560,420])
hold on, box on, grid on
for kS = 1:num_scans
    scatter(zClut{kS}(1,:),zClut{kS}(2,:),'.')
%     drawnow
end
axis equal
ax = gca;
ax.XTick = -Lx/2:lx:Lx/2;
ax.YTick = -Ly/2:ly:Ly/2;

figure('Position',[700,200,560,420])
imagesc(xc,yc,map_classic)
colorbar
axis xy
ax = gca;
ax.XTick = -Lx/2:lx:Lx/2;
ax.YTick = -Ly/2:ly:Ly/2;

figure('Position',[700,200,560,420])
imagesc(xc,yc,cm_obj.Map)
colorbar
axis xy
ax = gca;
ax.XTick = -Lx/2:lx:Lx/2;
ax.YTick = -Ly/2:ly:Ly/2;

% figure('Position',[700,200,560,420])
% imagesc(xc,yc,1./map_spatial)
% colorbar
% axis xy
% ax = gca;
% ax.XTick = -Lx/2:lx:Lx/2;
% ax.YTick = -Ly/2:ly:Ly/2;


