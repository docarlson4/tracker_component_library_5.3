function [xNew, SNew, lgclNew] = TwoPoint(obj, tCur, zCur, SRCur)
%TWO_POINT_INIT Two-point state and lower Cholesky decomposed covariance
%initialization for arbitrary spatial dimensions D
%
% INPUT
%   tCur - mmt time stamp (1 X N)
%   zCur - mmt vectors (D X N)
%   SRCur - mmt covariances (D X D X N)
%
% OUTPUT
%   xNew - state vectors
%   SNew - state lower Cholesky decomposed covariance
%
% USAGE
%   [xNew, SNew] = TwoPoint(obj, tCur, zCur, SRCur)

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2024-03-23 20:45

%%

persistent kBuf zBuf SBuf tBuf

vLims = [obj.VelMin, obj.VelMax];
xDim = obj.MotionModelObject.StateDim;
[zDim, numMeas] = size(zCur);

if isempty(kBuf)
    kBuf = 0;
    zBuf = cell(zDim,0);
    SBuf = cell(zDim,0);
    tBuf = cell(zDim,0);
end

xNew = zeros(xDim,0);
SNew = zeros(xDim,xDim,0);
lgclNew = zeros(1,0);

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
        t1 = tBuf{1}(1,:); N1 = length(t1);
        t2 = tBuf{2}(1,:); N2 = length(t2);
        del_t(1,:,:) = t2' - t1;

        % Mmts: D X N2 X N1, D - spatial dimension
        z1 = permute(repmat(zBuf{1},[1,1,N2]),[1,3,2]);
        z2 = repmat(zBuf{2},[1,1,N1]);
        del_z = z2 - z1;

        % Velocity and speed
        v = bsxfun(@rdivide,del_z,repmat(del_t,[2,1,1]));
        vMag = sqrt(squeeze(sum(v.*v,1)));

        % Apply speed limits
        lgclVlim = (vLims(1) < vMag) & (vMag < vLims(2));
        [i2,i1] = find(lgclVlim);
        
        if length(unique(i2)) < length(i2)
            disp("length(unique(i2)) < length(i2)")
        end

        % State vector
        sp = z2(:,i2);
        sv = zeros(size(sp));
        for k = 1:length(i1)
            sv(:,k) = v(:,i2(k),i1(k));
        end
        if obj.MotionModelObject.Type == "NCA"
            xNew = [sp;sv;zeros(size(sp))];
        elseif obj.MotionModelObject.Type == "NCV"
            xNew = [sp;sv];
        else
            error("Unknown Motion Model: " + obj.MotionModelObject.Type)
        end

        % State covariance
        detS = zeros(1,length(i1));
        for k = 1:length(i1)
            T = del_t(1,i2(k),i1(k));
            S22 = SBuf{1}(:,:,i1(k))/T;
            if obj.MotionModelObject.Type == "NCA"
                acc = obj.MotionModelObject.MaxKinVar ...
                    / sqrt(2 + zDim);
                SNew(:,:,k) = blkdiag([
                    SBuf{2}(:,:,i2(k)), zeros(zDim,zDim);
                    SBuf{2}(:,:,i2(k))/T, S22], ...
                    acc * eye(obj.SpaceDim));
            elseif obj.MotionModelObject.Type == "NCV"
                SNew(:,:,k) = [
                    SBuf{2}(:,:,i2(k)), zeros(zDim,zDim);
                    SBuf{2}(:,:,i2(k))/T, S22];
            else
                error("Unknown Motion Model: " + obj.MotionModelObject.Type)
            end
            detS(k) = det(SNew(:,:,k));
        end

        % Remove duplicates in favor of small detS
%         u2 = unique(i2); numU = length(u2);
%         if numU < length(i2)
%             iKeep = zeros(1,numU);
%             for kU = 1:numU
%                 ig = find(u2(kU)==i2);
%                 [~, iMinDetS] = min(detS);
%                 iKeep(kU) = ig(iMinDetS);
%             end
%             xNew = xNew(:,iKeep);
%             SNew = SNew(:,:,iKeep);
%         end

        lgclNew = false(1,numMeas);
        lgclNew(i2) = true;

        % Reset zBuf for next mmt
        zBuf{1} = zBuf{2};
        SBuf{1} = SBuf{2};
        tBuf{1} = tBuf{2};
        kBuf = 1;
    end
end

end
