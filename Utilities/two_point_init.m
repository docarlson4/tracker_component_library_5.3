function [s, S] = two_point_init(zCur, SRCur, tCur, vLims, xDim)
%TWO_POINT_INIT Summary of this function goes here
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
% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-23 20:45

%%

persistent kBuf zBuf SBuf tBuf

if isempty(kBuf)
    kBuf = 0;
    zBuf = cell(2,0);
    SBuf = cell(2,0);
    tBuf = cell(2,0);
end

s = zeros(xDim,0);
S = zeros(xDim,xDim,0);

numMeas = size(zCur,2);
if isscalar(tCur)
    tCur = tCur*ones(1,numMeas);
end
if size(SRCur,3) == 1
    SRCur = repmat(SRCur, [1,1,numMeas]);
end
if numMeas > 0
    kBuf = kBuf+1;
    zBuf{kBuf} = zCur;
    SBuf{kBuf} = SRCur;
    tBuf{kBuf} = tCur;
    if kBuf == 2
        % Generate initial state
        t1 = tBuf{1}(1,:);
        x1 = zBuf{1}(1,:);
        y1 = zBuf{1}(2,:);

        t2 = tBuf{2}(1,:);
        x2 = zBuf{2}(1,:);
        y2 = zBuf{2}(2,:);

        del_t = t2' - t1;
        del_x = x2' - x1;
        del_y = y2' - y1;
        vx = del_x./del_t;
        vy = del_y./del_t;
        v = sqrt(vx.*vx + vy.*vy);

        lgcl_v = (vLims(1) < v) & (v < vLims(2));

        [i2,i1] = find(lgcl_v);
        
        % State vector
        sx = x2(i2);
        sy = y2(i2);
        svx = zeros(size(sx));
        svy = zeros(size(sy));
        for k = 1:length(i1)
            svx(k) = vx(i2(k),i1(k));
            svy(k) = vy(i2(k),i1(k));
        end
        s = [sx;sy;svx;svy];

        % State covariance
        for k = 1:length(i1)
            T = del_t(i2(k),i1(k));
            S22 = SBuf{1}(:,:,i1(k))/T;
            S(:,:,k) = [
                SBuf{2}(:,:,i2(k)), zeros(2,2);
                SBuf{2}(:,:,i2(k))/T, S22];
        end

        % Reset zBuf for next mmt
        zBuf{1} = zBuf{2};
        SBuf{1} = SBuf{2};
        tBuf{1} = tBuf{2};
        kBuf = 1;
    end
end

end
