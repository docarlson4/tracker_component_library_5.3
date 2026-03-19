clc
clear
close all

cd(fileparts(which(mfilename)))

dbstop if error

MyConstants

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Doug Carlson (doug.o.carlson@gmail.com), 2024-03-23 20:02

%%
Vmin = 25*mph;
Vmax = 600*mph;
vLims = [Vmin; Vmax];

%%
load("swarm_mmts.mat","zMeasCart","SRMeasCart", "tMeas")

%%
figure('Position', [80,50,950,700])
% clf
hold on, grid on, box on

clear two_point_init
num_scan = length(zMeasCart);
for kS = 1:num_scan

    tCur = tMeas{kS};
    zCur = zMeasCart{kS};
    SRCur = SRMeasCart{kS};

    if(~isempty(zCur))
        scatter(zCur(1,:)/km,zCur(2,:)/km,'o');
    end

    [x, S] = two_point_init(zCur, SRCur, tCur, vLims);

    if ~isempty(x)
        scatter(x(1,:)/km,x(2,:)/km,'k.','linewidth',1)
    end

    h1 = xlabel('\bfx km');
    h2 = ylabel('\bfy km');
    h3 = title('\bf\fontsize{14}True and Estimated Trajectories');
    set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
    axis square

    axis([-1 1 -1 1]*50)

    drawnow
end



