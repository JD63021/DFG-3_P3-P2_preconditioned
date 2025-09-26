function uvp = build_residual_convective_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U1, U2, U3, gamma, mu, rho, t, bcFlags, inletProfile, corner) %#ok<INUSD>
% Nonlinear convection-only residual for P3–P2:
%   V = ∫ (u·∇u_x) N_i,  W = ∫ (u·∇u_y) N_i
% Uses the same 1/10 scaling pattern as your original code.
% BC handling here: ONLY zero convective residual at velocity Dirichlet DOFs.
% No pressure contribution from convection → P = 0.

% ---------- sizes ----------
numVelNodes = length(nodeInfo.velocity.x); % P3 velocity nodes
Nxy         = numVelNodes;
numElements = size(elemInfo.velElements,1);
Npr         = max(elemInfo.presElements(:)); %#ok<NASGU>

% ---------- triplets ----------
I_v = zeros(10*numElements,1); Val_v = zeros(10*numElements,1);
I_w = zeros(10*numElements,1); Val_w = zeros(10*numElements,1);
off_v = 0; off_w = 0;

% ---------- reference shapes ----------
[NV,dNxiV,dNetaV,wV,~] = precomputeShapeFunctionsP3_Tri();
numGauss = length(wV);

rhoRe = rho*Re;

% ---------- element loop ----------
for e = 1:numElements
    Kvel = elemInfo.velElements(e,:);

    xV = nodeInfo.velocity.x(Kvel);
    yV = nodeInfo.velocity.y(Kvel);

    U1el = U1(Kvel);
    U2el = U2(Kvel);

    vLocal = zeros(10,1);
    wLocal = zeros(10,1);

    % precompute shapes and geometry per GP
    shapeV  = zeros(10,numGauss);
    dNxV    = zeros(10,numGauss);
    dNyV    = zeros(10,numGauss);
    detJall = zeros(numGauss,1);

    for gp = 1:numGauss
        shapeV(:,gp) = NV(:,gp);
        [dNx,dNy,detJ] = p3ShapeDerivativesAllNodes(xV,yV,dNxiV(:,gp),dNetaV(:,gp));
        dNxV(:,gp) = dNx;  dNyV(:,gp) = dNy;  detJall(gp) = detJ;
    end

    % flow values at GPs
    a1s = zeros(numGauss,1); a2s = zeros(numGauss,1);
    du_dx_s = zeros(numGauss,1); du_dy_s = zeros(numGauss,1);
    dv_dx_s = zeros(numGauss,1); dv_dy_s = zeros(numGauss,1);

    for gp = 1:numGauss
        Ni   = shapeV(:,gp);
        dNxg = dNxV(:,gp); dNyg = dNyV(:,gp);

        a1s(gp)     = Ni.'*U1el;
        a2s(gp)     = Ni.'*U2el;
        du_dx_s(gp) = dNxg.'*U1el;  du_dy_s(gp) = dNyg.'*U1el;
        dv_dx_s(gp) = dNxg.'*U2el;  dv_dy_s(gp) = dNyg.'*U2el;
    end

    % convection residual (keep inner kNode loop with 1/10 scaling)
    for iNode = 1:10
        sumV_i = 0; sumW_i = 0;

        for kNode = 1:10 %#ok<NASGU> % kept to preserve 1/10 scaling pattern
            accumV = 0; accumW = 0;

            for gp = 1:numGauss
                wt   = wV(gp);
                detJ = detJall(gp);

                Ni_i = shapeV(iNode,gp);
                a1   = a1s(gp); a2 = a2s(gp);
                du_dx = du_dx_s(gp); du_dy = du_dy_s(gp);
                dv_dx = dv_dx_s(gp); dv_dy = dv_dy_s(gp);

                conv_u = a1*du_dx + a2*du_dy;
                conv_v = a1*dv_dx + a2*dv_dy;

                accumV = accumV + rhoRe*(1/10)*Ni_i*conv_u*(wt*detJ);
                accumW = accumW + rhoRe*(1/10)*Ni_i*conv_v*(wt*detJ);
            end

            sumV_i = sumV_i + accumV;
            sumW_i = sumW_i + accumW;
        end

        vLocal(iNode) = vLocal(iNode) + sumV_i;
        wLocal(iNode) = wLocal(iNode) + sumW_i;
    end

    % store local into triplets
    idxv = off_v+(1:10);  I_v(idxv) = Kvel(:);  Val_v(idxv) = vLocal;  off_v = off_v + 10;
    idxw = off_w+(1:10);  I_w(idxw) = Kvel(:);  Val_w(idxw) = wLocal;  off_w = off_w + 10;
end

% build global convection residual
V = accumarray(I_v,Val_v,[Nxy,1]);
W = accumarray(I_w,Val_w,[Nxy,1]);
P = zeros(max(elemInfo.presElements(:)),1); % convection gives no pressure residual

% zero out velocity BC rows (do NOT inject u-g here)
% Impose inlet target g via the convective residual so that K_lin identity rows yield (U - g)
[inletNodes, sideNodes] = collect_bc_sets_for_zero(boundaryInfo,bcFlags);

% Inlet: subtract the profile on x-velocity; inlet v_g = 0
if ~isempty(inletNodes)
    yAll = nodeInfo.velocity.y;
    H    = max(yAll) - min(yAll);
    yy   = yAll(inletNodes);
    UinVals = arrayfun(@(yyi) inletProfile(t, yyi, H), yy);
    V(inletNodes) = V(inletNodes) - UinVals;   % makes total residual row behave like (U_x - g)
    W(inletNodes) = 0;                         % v_g = 0 so no offset needed
end

% Walls/sides: homogeneous Dirichlet → zero the residual rows
if ~isempty(sideNodes)
    V(sideNodes) = 0;
    W(sideNodes) = 0;
end

uvp = [V; W; P];

end

% ===== helpers =====
function [inletNodes, sideNodes] = collect_bc_sets_for_zero(boundaryInfo,bcFlags)
% inlet
if isfield(bcFlags,'inlet'), topFlag = bcFlags.inlet; else, topFlag = 999999; end
if isfield(boundaryInfo, ['flag_' num2str(topFlag)])
    inletNodes = boundaryInfo.(['flag_' num2str(topFlag)])(:);
else
    inletNodes = [];
end
% walls/sides
if isfield(bcFlags,'wall'), sides = bcFlags.wall; else, sides = []; end
sideNodes = [];
for i = 1:length(sides)
    fn = ['flag_' num2str(sides(i))];
    if isfield(boundaryInfo,fn)
        sideNodes = [sideNodes; boundaryInfo.(fn)(:)]; %#ok<AGROW>
    end
end
sideNodes = unique(sideNodes);
end

function [dNxAll,dNyAll,detJ] = p3ShapeDerivativesAllNodes(x,y,dNxi,dNeta)
dX_dxi  = sum(x.*dNxi);   dX_deta = sum(x.*dNeta);
dY_dxi  = sum(y.*dNxi);   dY_deta = sum(y.*dNeta);
J = [dX_dxi dY_dxi; dX_deta dY_deta];
detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
invJ = inv(J);
dNxAll = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
dNyAll = invJ(2,1)*dNxi + invJ(2,2)*dNeta;
end
