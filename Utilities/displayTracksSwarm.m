%% Plots
disp('Displaying contacts, true tracks, and tracks found.') 

figure('Position', [80,50,950,700])
% clf
hold on, grid on, box on

PDisp = 0.98;

trackList = AVLTree();
for curScan = 1:numSamples
    if ~isempty(ZCartTrue{curScan})
        xTrue = ZCartTrue{curScan};
        scatter(xTrue(2,:)/km,xTrue(1,:)/km,'bo','linewidth',1)
    end
    
    zCur = zMeasCart{curScan};
    if(~isempty(zCur))
        scatter(zCur(1,:)/km,zCur(2,:)/km,'ok');
    end
    
    stateCur = xStates{curScan,1};
    IDCur = xStates{curScan,2};
    rCur = rStates{curScan};
    
    numHyp = length(rCur);
    
    trackListNew = AVLTree();
    
    for curHyp = 1:numHyp
        curID = IDCur(curHyp);
        prevTrackLoc = trackList.find(curID);
        
        %If it is an existing track, then connect the points
        if(~isempty(prevTrackLoc))
            prevTrackLoc = prevTrackLoc.value;
            
            curTrackLoc = stateCur(1:2,curHyp);
            plot([prevTrackLoc(1);curTrackLoc(1)]/km,[prevTrackLoc(2);curTrackLoc(2)]/km,'-g','linewidth',2);
            trackListNew.insert(KeyVal(curID,curTrackLoc));
        elseif(rCur(curHyp)>PDisp)
            %We only start the tentative track when it is above the
            %detection threshold for track existence.
            curTrackLoc = stateCur(1:2,curHyp);
            
            if(~isempty(prevTrackLoc))
                prevTrackLoc = prevTrackLoc.value;
                plot([prevTrackLoc(1);curTrackLoc(1)]/km,[prevTrackLoc(2);curTrackLoc(2)]/km,'-g','linewidth',2);
            else
                scatter(curTrackLoc(1)/km,curTrackLoc(2)/km,'.g')
            end
            
            trackListNew.insert(KeyVal(curID,curTrackLoc));
        end
    end
    
    trackList = trackListNew;


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

