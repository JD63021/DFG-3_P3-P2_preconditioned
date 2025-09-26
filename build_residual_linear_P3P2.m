function uvp = build_residual_linear_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U1, U2, U3, gamma, mu, rho, t, bcFlags, inletProfile, corner)
% Linear (u-independent in coefficients) residual for P3–P2:
%  - viscosity + grad-div (velocity)
%  - pressure gradient (momentum) and divergence (continuity)
%  - Dirichlet BCs enforced as residual: u - g = 0 on velocity BCs
%  - pressure pin enforced as residual: p(corner) - 0 = 0
%
% Convection NOT included here.

% ---------- sizes ----------
numVelNodes = length(nodeInfo.velocity.x); % P3 velocity nodes
Nxy         = numVelNodes;
numElements = size(elemInfo.velElements,1);
Npr         = max(elemInfo.presElements(:));

% ---------- triplets for element -> global residual ----------
I_v = zeros(10*numElements,1); Val_v = zeros(10*numElements,1);
I_w = zeros(10*numElements,1); Val_w = zeros(10*numElements,1);
I_p = zeros( 6*numElements,1); Val_p = zeros( 6*numElements,1);
off_v=0; off_w=0; off_p=0;

% ---------- inlet geometry (for profile) ----------
yMin = min(nodeInfo.velocity.y);
yMax = max(nodeInfo.velocity.y);
H    = yMax - yMin;

% ---------- reference shapes ----------
[NV,dNxiV,dNetaV,wV,~] = precomputeShapeFunctionsP3_Tri();
[NP,dNxiP,dNetaP,~,~]  = precomputeShapeFunctionsP2_Tri();
numGauss = length(wV);

% ---------- element loop ----------
for e=1:numElements
    Kvel = elemInfo.velElements(e,:);   % 10 P3 nodes
    Kpr  = elemInfo.presElements(e,:);  %  6 P2 nodes

    xV = nodeInfo.velocity.x(Kvel);  yV = nodeInfo.velocity.y(Kvel);

    % local solution
    U1el = U1(Kvel);  U2el = U2(Kvel);  Pel = U3(Kpr);

    % local accumulators
    vLocal = zeros(10,1);
    wLocal = zeros(10,1);
    pLocal = zeros( 6,1);

    % precompute shapes & detJ for all gps
    shapeV  = zeros(10,numGauss);
    dNxV    = zeros(10,numGauss);
    dNyV    = zeros(10,numGauss);
    detJall = zeros(numGauss,1);
    shapeP  = zeros(6,numGauss);

    for gp=1:numGauss
        shapeV(:,gp) = NV(:,gp);
        [dNx,dNy,detJ] = p3ShapeDerivativesAllNodes( ...
                            xV,yV,dNxiV(:,gp),dNetaV(:,gp));
        dNxV(:,gp) = dNx;  dNyV(:,gp) = dNy;  detJall(gp) = detJ;
        if gp <= size(NP,2)
            shapeP(:,gp) = NP(:,gp);
        end
    end

    % -------- velocity residual: viscosity + grad-div + pressure gradient
    for iNode=1:10
        sumV_i = 0; sumW_i = 0;

        for kNode=1:10
            accumV = 0; accumW = 0;

            for gp=1:numGauss
                wt   = wV(gp);
                detJ = detJall(gp);

                dNx_i = dNxV(iNode,gp);   dNy_i = dNyV(iNode,gp);
                dNx_k = dNxV(kNode,gp);   dNy_k = dNyV(kNode,gp);

                Uk1 = U1el(kNode);        Uk2 = U2el(kNode);

                % viscosity + grad-div (component form)
                vx_xx = mu*(2*dNx_i*dNx_k + dNy_i*dNy_k) + gamma*(dNx_i*dNx_k);
                vx_xy = mu*(dNy_i*dNx_k)                 + gamma*(dNx_i*dNy_k);
                vy_yx = mu*(dNx_i*dNy_k)                 + gamma*(dNy_i*dNx_k);
                vy_yy = mu*(dNx_i*dNx_k + 2*dNy_i*dNy_k) + gamma*(dNy_i*dNy_k);

                accumV = accumV + ( Uk1*vx_xx + Uk2*vx_xy )*(wt*detJ);
                accumW = accumW + ( Uk1*vy_yx + Uk2*vy_yy )*(wt*detJ);
            end
            sumV_i = sumV_i + accumV;
            sumW_i = sumW_i + accumW;
        end

        % pressure gradient terms
        sumVp = 0; sumWp = 0;
        for pNode=1:6
            pVal = Pel(pNode);
            accumVp = 0; accumWp = 0;

            for gp=1:numGauss
                wt   = wV(gp);
                detJ = detJall(gp);

                dNx_i = dNxV(iNode,gp);  dNy_i = dNyV(iNode,gp);
                Np    = shapeP(pNode,gp);

                accumVp = accumVp - dNx_i * pVal * Np * (wt*detJ);
                accumWp = accumWp - dNy_i * pVal * Np * (wt*detJ);
            end
            sumVp = sumVp + accumVp;
            sumWp = sumWp + accumWp;
        end

        vLocal(iNode) = vLocal(iNode) + sumV_i + sumVp;
        wLocal(iNode) = wLocal(iNode) + sumW_i + sumWp;
    end

    % -------- continuity residual: -div(u) --------
    for pNode=1:6
        sumP = 0;
        for kNode=1:10
            Uk1 = U1el(kNode);  Uk2 = U2el(kNode);
            accumP = 0;

            for gp=1:numGauss
                wt   = wV(gp);
                detJ = detJall(gp);

                Np   = shapeP(pNode,gp);
                dNxk = dNxV(kNode,gp);  dNyk = dNyV(kNode,gp);

                accumP = accumP - (Uk1*dNxk + Uk2*dNyk) * Np * (wt*detJ);
            end
            sumP = sumP + accumP;
        end
        pLocal(pNode) = pLocal(pNode) + sumP;
    end

    % -------- store local into triplets --------
    idxv = off_v+(1:10);  I_v(idxv)=Kvel(:);  Val_v(idxv)=vLocal;  off_v=off_v+10;
    idxw = off_w+(1:10);  I_w(idxw)=Kvel(:);  Val_w(idxw)=wLocal;  off_w=off_w+10;
    idxp = off_p+(1:6);   I_p(idxp)=Kpr(:);   Val_p(idxp)=pLocal;  off_p=off_p+6;
end

% ---------- build global residual pieces ----------
V = accumarray(I_v,Val_v,[Nxy,1]);
W = accumarray(I_w,Val_w,[Nxy,1]);
P = accumarray(I_p,Val_p,[Npr,1]);

% ---------- apply Dirichlet BCs as residual: u - g = 0 ----------
[inletNodes, sideNodes, UinVals] = collect_bc_sets(nodeInfo,boundaryInfo,bcFlags,inletProfile,t);

% % inlet (u=profile, v=0 implied by your scheme)
% if ~isempty(inletNodes)
%     V(inletNodes) = U1(inletNodes) - UinVals;
%     W(inletNodes) = U2(inletNodes);
% end

% inlet handled in the convective residual to avoid rebuilding a BC vector each step.
% Do NOT touch inlet rows here (leave V/W as assembled):
% (no assignment to V(inletNodes) or W(inletNodes))


% walls/sides (u=0, v=0)
if ~isempty(sideNodes)
    V(sideNodes) = U1(sideNodes);
    W(sideNodes) = U2(sideNodes);
end

% ---------- pressure pin as residual ----------
P(corner) = U3(corner);

% ---------- pack ----------
uvp = [V; W; P];

end

% ===== helpers =====
function [inletNodes, sideNodes, UinVals] = collect_bc_sets(nodeInfo,boundaryInfo,bcFlags,inletProfile,t)
% inlet set
if isfield(bcFlags,'inlet'), topFlag = bcFlags.inlet; else, topFlag = 999999; end
if isfield(boundaryInfo, ['flag_' num2str(topFlag)])
    inletNodes = boundaryInfo.(['flag_' num2str(topFlag)])(:);
else
    inletNodes = [];
end
yMin = min(nodeInfo.velocity.y); yMax = max(nodeInfo.velocity.y); H = yMax - yMin;
if ~isempty(inletNodes)
    yy = nodeInfo.velocity.y(inletNodes);
    UinVals = arrayfun(@(yyi) inletProfile(t,yyi,H), yy);
else
    UinVals = [];
end

% wall/side set
if isfield(bcFlags,'wall'), sides = bcFlags.wall; else, sides = []; end
sideNodes = [];
for i=1:length(sides)
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
