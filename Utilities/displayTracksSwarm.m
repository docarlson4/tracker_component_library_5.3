%% Plots
disp('Displaying contacts, true tracks, and tracks found.') 

figure('Position', [80,50,950,700])
% clf
hold on, grid on, box on

ttl(1) = "\bf\fontsize{14}True and Estimated Trajectories";
ttl(2) = "\bf\fontsize{12}Moton Model: " + motion_model ...
    + ", State Initialization: " + init_method;

PDisp = 0.98;

trackList = AVLTree();
for curScan = 1:numSamples

    % Plot truth data
    if ~isempty(ZCartTrue{curScan})
        xTrue = ZCartTrue{curScan};
        scatter(xTrue(2,:)/km,xTrue(1,:)/km,'bo','linewidth',1)
    end
    
    % Plot measurements
    zCur = zMeasCart{curScan};
    if(~isempty(zCur))
        scatter(zCur(1,:)/km,zCur(2,:)/km,'ok');
    end
    
    xCur = state_st(curScan).x;
    SCur = state_st(curScan).S;
    IDCur = state_st(curScan).ID;
    rCur = state_st(curScan).r;
    
    numHyp = length(rCur);

    trackListNew = AVLTree();
    
    for curHyp = 1:numHyp
        curID = IDCur(curHyp);
        prevTrackLoc = trackList.find(curID);
        
        %If it is an existing track, then connect the points
        if(~isempty(prevTrackLoc))
            prevTrackLoc = prevTrackLoc.value;
            
            curTrackLoc = xCur(:,curHyp);
            plot([prevTrackLoc(1);curTrackLoc(1)]/km,[prevTrackLoc(2);curTrackLoc(2)]/km,'-g','linewidth',2);
            trackListNew.insert(KeyVal(curID,curTrackLoc));
        elseif(rCur(curHyp)>PDisp)
            %We only start the tentative track when it is above the
            %detection threshold for track existence.
            curTrackLoc = xCur(:,curHyp);
            
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
    h3 = title(ttl);
    set(gca,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h1,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h2,'FontSize',14,'FontWeight','bold','FontName','Times')
    set(h3,'FontSize',14,'FontWeight','bold','FontName','Times')
    axis square

    axis([-1 1 -1 1]*50)

    drawnow
end

%% Plot Tracks 

% Filter minimum existence probability: PIsRealInit
rCell = {track_st.r}';
lgcl_r = cellfun(@(x) any(x > 0.5), rCell, 'UniformOutput',false);
lgcl_r = cell2mat(lgcl_r);
track_st = track_st(lgcl_r);
% error("double check Filter minimum existence probability")

% Filters minimum number of hits: num_hits
min_hits = 15;
num_hits = [track_st.num_hits]';
lgcl_hits = num_hits >= min_hits;
track_st = track_st(lgcl_hits);

% Plot remaining tracks
num_trx = length(track_st);
clrs = lines(num_trx);
line_style = ["p-","*-","s-","v-","^-","d-"];
lgcl_hp = false(1,num_trx);
for kT = num_trx:-1:1
    x  = track_st(kT).x(1,:);
    y  = track_st(kT).x(2,:);
    vx = track_st(kT).x(3,:);
    vy = track_st(kT).x(4,:);
    if any(lgcl_r)
        kL = floor((kT-1)/7)+1;
        hp(kT) = plot(x/km,y/km,line_style(kL),'Color',clrs(kT,:), ...
            LineWidth=2, DisplayName="track "+num2str(kT));
        lgcl_hp(kT) = true;
    end
end
legend(hp(lgcl_hp))



