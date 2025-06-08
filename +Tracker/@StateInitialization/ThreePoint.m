function [xNew, SNew, lgclNew] = ThreePoint(obj, tCur, zCur, SRCur)
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

% Douglas O. Carlson, Ph.D. (doug.o.carlson@gmail.com), 2025-06-07

%%
vLims = [obj.VelMin, obj.VelMax];
accMax = obj.AccMax;
geeMax = obj.GeeMax;
cosPhiMax = cos(obj.PhiMax);
xDim = obj.MotionModelObject.StateDim;
[zDim, numMeas] = size(zCur);
Z00 = zeros(zDim,zDim);
zD = zeros(3-zDim,1);

lgclNew = false(1,numMeas);

xNew = zeros(xDim,0);
SNew = zeros(xDim,xDim,0);

if isscalar(tCur)
    tCur = tCur*ones(1,numMeas);
end
if size(SRCur,3) == 1
    SRCur = repmat(SRCur, [1,1,numMeas]);
end
if numMeas > 0
    obj.kBuf = obj.kBuf+1;
    obj.zBuf{obj.kBuf} = zCur;
    obj.SBuf{obj.kBuf} = SRCur;
    obj.tBuf{obj.kBuf} = tCur;
    if obj.kBuf == 3
        % Generate initial state
        t1 = obj.tBuf{1}(1,:); N1 = length(t1);
        t2 = obj.tBuf{2}(1,:); N2 = length(t2);
        t3 = obj.tBuf{3}(1,:); N3 = length(t3);
        del_t21(1,:,:) = t2' - t1;
        del_t32(1,:,:) = t3' - t2;

        % Mmts: D X N3 X N1, D - spatial dimension
        z1 = permute(repmat(obj.zBuf{1},[1,1,N2]),[1,3,2]);
        z2 = repmat(obj.zBuf{2},[1,1,N1]);
        del_z21 = z2 - z1;
        z2 = permute(repmat(obj.zBuf{2},[1,1,N3]),[1,3,2]);
        z3 = repmat(obj.zBuf{3},[1,1,N2]);
        del_z32 = z3 - z2;

        % Velocity and speed
        v21 = bsxfun(@rdivide,del_z21,repmat(del_t21,[2,1,1]));
        vMag21 = sqrt(sum(v21.*v21,1));
        sqvMag21 = mySqueeze(vMag21);
        vHat21 = bsxfun(@rdivide,v21,vMag21);
        v32 = bsxfun(@rdivide,del_z32,repmat(del_t32,[2,1,1]));
        vMag32 = sqrt(sum(v32.*v32,1));
        sqvMag32 = mySqueeze(vMag32);
        vHat32 = bsxfun(@rdivide,v32,vMag32);

        % Apply speed limits
        lgcl_vMag21 = sparse((vLims(1) < sqvMag21) & (sqvMag21 < vLims(2)));
        lgcl_vMag32 = sparse((vLims(1) < sqvMag32) & (sqvMag32 < vLims(2)));

        % Apply max-acceleration filter
        % aMag = sqrt(mySqueeze(sum((v32-v21).*(v32-v21),1)));
        % lgcl_aMax = (aMag < accMax)

        idx_triplets = FindConnectedPaths(lgcl_vMag21, lgcl_vMag32);
        if ~isempty(idx_triplets)
            lgclNew(idx_triplets(:,3)) = true;
        end

        sp = []; sv = []; sa = [];
        for k = 1:size(idx_triplets, 1)
            % State vector
            % - position
            % - velocity
            v = v32(:, idx_triplets(k,3), idx_triplets(k,2));
            % - acceleration
            a = v - v21(:, idx_triplets(k,2), idx_triplets(k,1));
            aDen = del_t32(1, idx_triplets(k,3), idx_triplets(k,2)) ...
                + del_t21(1, idx_triplets(k,2), idx_triplets(k,1));
            a = a./aDen;
            % - angle limit: c = cos(phi) == vH32.vH21
            c = dot(vHat32(:,idx_triplets(k,3), idx_triplets(k,2)), ...
                vHat21(:,idx_triplets(k,2), idx_triplets(k,1)));
            % - g-force limit
            vM = (v + v21(:, idx_triplets(k,2), idx_triplets(k,1)))/2;
            ac = vecnorm(cross([vM;zD], [a;zD]))./vecnorm(vM);
            g = ac/9.81;
            if vecnorm(a) > accMax || abs(c) < cosPhiMax || g > geeMax
                lgclNew(idx_triplets(k,3)) = false;
                continue;
            end
            sp = cat( 2, sp, z3(:, idx_triplets(k,3)) );
            sv = cat( 2, sv, v );
            sa = cat( 2, sa, a );

            % State covariance
            T21 = del_t21(1, idx_triplets(k,2), idx_triplets(k,1));
            T32 = del_t32(1, idx_triplets(k,3), idx_triplets(k,2));
            T31 = T32 + T21;

            S11 = obj.SBuf{3}(:, :, idx_triplets(k,3));
            S21 = S11/T32;
            S31 = S11/(T32*T31);
            S22 = obj.SBuf{2}(:, :, idx_triplets(k,2))/T32;
            S32 = S22/T21;
            S33 = obj.SBuf{1}(:, :, idx_triplets(k,1))/(T31*T21);

            SNew(:,:,k) = [
                S11, Z00, Z00
                S21, S31, Z00;
                S31, S32, S33];
        end
        xNew = [sp;sv;sa];

        % Reset zBuf for next mmt
        obj.zBuf{1} = obj.zBuf{2};
        obj.SBuf{1} = obj.SBuf{2};
        obj.tBuf{1} = obj.tBuf{2};
        obj.zBuf{2} = obj.zBuf{3};
        obj.SBuf{2} = obj.SBuf{3};
        obj.tBuf{2} = obj.tBuf{3};
        obj.kBuf = 2;
    end
end

end
        %         % Angle constraint
        %         cos_phi0 = cosd(90);
        %         for k2 = 1:N2
        %             del_z32 = z3(:,i3) - z2(:,k2);
        %             mag_z32 = vecnorm(del_z32);
        %             del_z21 = z2(:,k2) - z1(:,i1);
        %             mag_z21 = vecnorm(del_z21);
        %             cos_phi = (dot(del_z32,del_z21)./(mag_z32.*mag_z21 + eps));
        %             lgcl_phi = cos_phi > cos_phi0;
        %         end

