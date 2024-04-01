function [track_st] = get_tracks(curScan, state_st, track_st)
%GET_TRACKS Summary of this function goes here
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

% Developed in Matlab 9.12.0.2529717 (R2022a) Update 8 on PCWIN64.
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-31 20:35

%%

xCur = state_st.x;
SCur = state_st.S;
IDCur = state_st.ID;
rCur = state_st.r;

IDtrx = [track_st.ID];
if isempty(IDtrx) && ~isempty(IDCur)
    numHyp = length(IDCur);
    track_st = repmat(track_st,[numHyp,1]);
    for kH = 1:numHyp
        track_st(kH).x = xCur(:,kH);
        track_st(kH).S = SCur(:,:,kH);
        track_st(kH).ID = IDCur(kH);
        track_st(kH).r = rCur(kH);
        track_st(kH).scan_num = curScan;
        track_st(kH).num_hits = 1;
    end
end

if ~isempty(IDtrx) && ~isempty(IDCur)
    [IDOld, idx_old, idx_cur] = intersect(IDtrx, IDCur);
    for kO = 1:length(IDOld)
        idO = idx_old(kO);
        track_st(idO).x = [track_st(idO).x, xCur(:,idx_cur(kO))];
        track_st(idO).S = cat(3, track_st(idO).S, SCur(:,:,idx_cur(kO)));
        % track_st(kO).ID = [track_st(idO)).ID, IDCur(kH)];
        track_st(idO).r = [track_st(idO).r, rCur(idx_cur(kO))];
        track_st(idO).scan_num = [track_st(idO).scan_num, curScan];
        track_st(idO).num_hits = track_st(idO).num_hits + 1;
    end
    [~, idx_new] = setdiff(IDCur, IDOld);
    numTrx = length(track_st);
    for kN = idx_new'
        track_st(kN+numTrx).x = xCur(:,kN);
        track_st(kN+numTrx).S = SCur(:,:,kN);
        track_st(kN+numTrx).ID = IDCur(kN);
        track_st(kN+numTrx).r = rCur(kN);
        track_st(kN+numTrx).scan_num = curScan;
        track_st(kN+numTrx).num_hits = 1;
    end
end



end
