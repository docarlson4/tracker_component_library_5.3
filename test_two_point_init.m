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

num_scan = length(zMeasCart);
zBuf = cell(2,1);
SBuf = cell(2,1);
tBuf = cell(2,1);

kBuf = 0;
for kS = 1:num_scan

    zCur = zMeasCart{kS};
    if(~isempty(zCur))
        scatter(zCur(1,:)/km,zCur(2,:)/km,'o');
    end

    num_mmt = size(zMeasCart{kS},2);
    if num_mmt > 0
        kBuf = kBuf+1;
        zBuf{kBuf} = zMeasCart{kS};
        SBuf{kBuf} = SRMeasCart{kS};
        tBuf{kBuf} = tMeas{kS};
        if kBuf == 2
            % Generate initial state
            [x, S] = two_point_init(zBuf, SBuf, tBuf, vLims);

            if ~isempty(x)
                scatter(x(1,:)/km,x(2,:)/km,'k.','linewidth',1)
            end

            % Reset zBuf for next mmt
            zBuf{1} = zBuf{2};
            SBuf{1} = SBuf{2};
            tBuf{1} = tBuf{2};
            kBuf = 1;
        end
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



