function uvp = build_residual_convective_P3P2_fast( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U1, U2, U3, gamma, mu, rho, t, bcFlags, inletProfile, corner) %#ok<INUSD,INUSL>
% Convective residual only (P3–P2), vectorized, no Gauss loops.
% Returns [V; W; zeros(Npr,1)] with inlet BC injected here (u - g).
%
% NOTE: Linear (viscous/pressure) parts are excluded on purpose.

% sizes
Nxy = length(nodeInfo.velocity.x);
Npr = max(elemInfo.presElements(:));
nEl = size(elemInfo.velElements,1);

% convective scaling (match your P3–P2 residual scaling)
beta = rho * Re * (1/1);

% geometry cache (built once, reused)
cache = get_cache_p3p2(nodeInfo, elemInfo);   % fields: NV, dNx{e}, dNy{e}, wdet{e}

V = zeros(Nxy,1);
W = zeros(Nxy,1);

for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);  % 1×10
    U1el = U1(Kvel);                   % 10×1
    U2el = U2(Kvel);

    NV   = cache.NV;                   % 10×nGp
    dNx  = cache.dNx{e};               % 10×nGp
    dNy  = cache.dNy{e};               % 10×nGp
    wg   = cache.wdet{e}(:);           % nGp×1

    % advecting field & grads at GP (vectorized)
    a1     = dgemvT(NV, U1el);         % nGp×1
    a2     = dgemvT(NV, U2el);
    du_dx  = dgemvT(dNx, U1el);
    du_dy  = dgemvT(dNy, U1el);
    dv_dx  = dgemvT(dNx, U2el);
    dv_dy  = dgemvT(dNy, U2el);

    conv_u = a1 .* du_dx + a2 .* du_dy;   % nGp×1
    conv_v = a1 .* dv_dx + a2 .* dv_dy;   % nGp×1

    % vLocal = ∫ Ni * conv_u * wdet
    vLocal = beta * ( NV * (wg .* conv_u) );   % 10×1
    wLocal = beta * ( NV * (wg .* conv_v) );   % 10×1

    V(Kvel) = V(Kvel) + vLocal;
    W(Kvel) = W(Kvel) + wLocal;
end

% --- BC handling in convective residual (as discussed) -------------------
[inletNodes, sideNodes] = get_bc_sets(boundaryInfo, bcFlags);

% inlet: impose (U - g) via residual offset on x-comp; v_g = 0
if ~isempty(inletNodes)
    yAll = nodeInfo.velocity.y;
    H    = max(yAll) - min(yAll);
    yy   = yAll(inletNodes);
    Uin  = arrayfun(@(yyi) inletProfile(t, yyi, H), yy);
    V(inletNodes) = V(inletNodes) - Uin(:);
    W(inletNodes) = 0;
end

% walls/sides: homogeneous Dirichlet → zero convective residual rows
if ~isempty(sideNodes)
    V(sideNodes) = 0;
    W(sideNodes) = 0;
end

P = zeros(Npr,1);   % no convective pressure part
uvp = [V; W; P];
end

% ===================== helpers (local) =====================

function cache = get_cache_p3p2(nodeInfo, elemInfo)
% Persistent geometry cache keyed by mesh sizes; rebuilt on mesh change.
persistent C Nxy_s nEl_s
Nxy  = length(nodeInfo.velocity.x);
nEl  = size(elemInfo.velElements,1);

if isempty(C) || ~isequal(Nxy, Nxy_s) || ~isequal(nEl, nEl_s)
    [NV,dNxiV,dNetaV,wV] = precomputeShapeFunctionsP3_Tri(); % 10×nGp
    nGp = numel(wV);
    C          = struct;
    C.NV       = NV;
    C.dNx      = cell(nEl,1);
    C.dNy      = cell(nEl,1);
    C.wdet     = cell(nEl,1);
    for e=1:nEl
        Kvel = elemInfo.velElements(e,:);
        xV = nodeInfo.velocity.x(Kvel);  yV = nodeInfo.velocity.y(Kvel);
        dNx = zeros(10,nGp); dNy = zeros(10,nGp); wdet = zeros(1,nGp);
        for gp=1:nGp
            [dNx(:,gp), dNy(:,gp), detJ] = p3ShapeDerivativesAllNodes( ...
                xV, yV, dNxiV(:,gp), dNetaV(:,gp));
            wdet(gp) = wV(gp) * detJ;
        end
        C.dNx{e}  = dNx;
        C.dNy{e}  = dNy;
        C.wdet{e} = wdet;
    end
    Nxy_s = Nxy; nEl_s = nEl;
end
cache = C;
end

function y = dgemvT(A, x)
% y = A.' * x  (A: m×n, x: m×1) → y: n×1
y = (A.' * x);
end

function [inletNodes, sideNodes] = get_bc_sets(boundaryInfo, bcFlags)
% Light helper; assumes flags exist when set.
if isfield(bcFlags,'inlet') && isfield(boundaryInfo, ['flag_' num2str(bcFlags.inlet)])
    inletNodes = boundaryInfo.(['flag_' num2str(bcFlags.inlet)])(:);
else
    inletNodes = [];
end
sideNodes = [];
if isfield(bcFlags,'wall')
    for f = reshape(bcFlags.wall,1,[])
        nm = ['flag_' num2str(f)];
        if isfield(boundaryInfo,nm)
            sideNodes = [sideNodes; boundaryInfo.(nm)(:)]; %#ok<AGROW>
        end
    end
    sideNodes = unique(sideNodes);
end
end