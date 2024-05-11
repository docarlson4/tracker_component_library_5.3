function [xNew, SNew] = ThreePoint(obj, tCur, zCur, SRCur)
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
xDim = obj.StateDim;
[zDim, numMeas] = size(zCur);

if isempty(kBuf)
    kBuf = 0;
    zBuf = cell(3,0);
    SBuf = cell(3,0);
    tBuf = cell(3,0);
end

xNew = zeros(xDim,0);
SNew = zeros(xDim,xDim,0);

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
    if kBuf == 3
        % Generate initial state
        t1 = tBuf{1}(1,:); N1 = length(t1);
        t2 = tBuf{2}(1,:); N2 = length(t2);
        t3 = tBuf{3}(1,:); N3 = length(t3);
        del_t21(1,:,:) = t2' - t1;
        del_t32(1,:,:) = t3' - t2;
        del_t31(1,:,:) = t3' - t1;

        % Mmts: D X N3 X N1, D - spatial dimension
        z1 = permute(repmat(zBuf{1},[1,1,N3]),[1,3,2]);
        z2 = zBuf{2};
        z3 = repmat(zBuf{3},[1,1,N1]);
        del_z31 = z3 - z1;

        % Velocity and speed
        v = bsxfun(@rdivide,del_z31,repmat(del_t31,[2,1,1]));
        vMag = sqrt(squeeze(sum(v.*v,1)));

        % Apply speed limits
        lgcl_vMag = (vLims(1) < vMag) & (vMag < vLims(2));
        [i3,i1] = find(lgcl_vMag);
        % num31 = length(i1);

        % Angle constraint
        cos_phi0 = cosd(90);
        for k2 = 1:N2
            del_z32 = z3(:,i3) - z2(:,k2);
            mag_z32 = vecnorm(del_z32);
            del_z21 = z2(:,k2) - z1(:,i1);
            mag_z21 = vecnorm(del_z21);
            cos_phi = (dot(del_z32,del_z21)./(mag_z32.*mag_z21 + eps));
            lgcl_phi = cos_phi > cos_phi0
        end
        
        % State vector
        sp = z3(:,i3);
        sv = zeros(size(sp));
        for k = 1:length(i1)
            sv(:,k) = v(:,i3(k),i1(k));
        end
        xNew = [sp;sv];

        % State covariance
        for k = 1:length(i1)
            T = del_t31(1,i3(k),i1(k));
            S22 = SBuf{1}(:,:,i1(k))/T;
            SNew(:,:,k) = [
                SBuf{3}(:,:,i3(k)), zeros(zDim,zDim);
                SBuf{3}(:,:,i3(k))/T, S22];
        end

        % Reset zBuf for next mmt
        zBuf{1} = zBuf{2};
        SBuf{1} = SBuf{2};
        tBuf{1} = tBuf{2};
        zBuf{2} = zBuf{3};
        SBuf{2} = SBuf{3};
        tBuf{2} = tBuf{3};
        kBuf = 1;
    end
end

end
