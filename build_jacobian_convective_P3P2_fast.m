function [Kc, Jx1, Jx2, Jy1, Jy2, Jx3, Jy3, Jc1, Jc2] = ...
    build_jacobian_convective_P3P2_fast( ...
      nodeInfo, elemInfo, boundaryInfo, ...
      D, Re, o, U1, U2, U3, gamma, mu, rho, corner) %#ok<INUSD,INUSL>
% Convective Jacobian only (P3–P2), vectorized, no Gauss loops.
% Outputs VV blocks (Jx1,Jx2,Jy1,Jy2); VP/PV/PP are zeros here.
% BC: zero Dirichlet velocity rows ONLY (no diag=1), no pressure pin.

Nxy = length(nodeInfo.velocity.x);
Npr = max(elemInfo.presElements(:));
nEl = size(elemInfo.velElements,1);
nV  = 10;

beta = rho * Re * (1/1);    % match residual scaling

% geometry cache
cache = get_cache_p3p2(nodeInfo, elemInfo);

% persistent triplet pattern for VV blocks (I,J and element ranges)
pat = get_conv_pattern(elemInfo, Nxy);

nn = numel(pat.I_VV);          % entries per VV block set
V_jx1 = zeros(nn,1);
V_jx2 = zeros(nn,1);
V_jy1 = zeros(nn,1);
V_jy2 = zeros(nn,1);

blk = nV*nV;

for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);  %#ok<NASGU> (indices implicit in pattern)

    NV   = cache.NV;          % 10×nGp
    dNx  = cache.dNx{e};      % 10×nGp
    dNy  = cache.dNy{e};      % 10×nGp
    wg   = cache.wdet{e}(:);  % nGp×1

    U1el = U1(elemInfo.velElements(e,:));  % 10×1
    U2el = U2(elemInfo.velElements(e,:));

    % advecting field & grads at GP
    a1     = dgemvT(NV, U1el);   % nGp×1
    a2     = dgemvT(NV, U2el);
    a1d1   = dgemvT(dNx, U1el);
    a1d2   = dgemvT(dNy, U1el);
    a2d1   = dgemvT(dNx, U2el);
    a2d2   = dgemvT(dNy, U2el);

    % helper: scale NV columns by (wg .* scalar_gp) then * NV' or dN' (GEMM)
    % A = NV .* (ones(nV,1) * (wg.*s(:)).')  → 10×nGp
    A_xx = (NV .* (ones(nV,1) * (wg.*a1d1).'));   % for Ni*Nk*a1d1
    A_xy = (NV .* (ones(nV,1) * (wg.*a1d2).'));
    A_yx = (NV .* (ones(nV,1) * (wg.*a2d1).'));
    A_yy = (NV .* (ones(nV,1) * (wg.*a2d2).'));

    Tdx  = (NV .* (ones(nV,1) * (wg.*a1).')) * dNx.';   % Ni*(a1 dNx_k)
    Tdy  = (NV .* (ones(nV,1) * (wg.*a2).')) * dNy.';   % Ni*(a2 dNy_k)

    % VV convective sub-blocks (Newton linearization)
    Kx1e = beta * ( A_xx * NV.' + Tdx + Tdy );   % ∫ Ni*(Nk*a1d1 + a·∇Nk)
    Kx2e = beta * ( A_xy * NV.' );               % ∫ Ni*(Nk*a1d2)
    Ky1e = beta * ( A_yx * NV.' );               % ∫ Ni*(Nk*a2d1)
    Ky2e = beta * ( A_yy * NV.' + Tdx + Tdy );   % ∫ Ni*(Nk*a2d2 + a·∇Nk)

    % scatter into persistent pattern slots
   % new (force double indices)
start = double(pat.offVV(e)) + 1;
rng   = start : (start + blk - 1);

    V_jx1(rng) = Kx1e(:);
    V_jx2(rng) = Kx2e(:);
    V_jy1(rng) = Ky1e(:);
    V_jy2(rng) = Ky2e(:);
end

% assemble sparse VV blocks
I = double(pat.I_VV); J = double(pat.J_VV);
Jx1 = sparse(I, J, V_jx1, Nxy, Nxy);
Jx2 = sparse(I, J, V_jx2, Nxy, Nxy);
Jy1 = sparse(I, J, V_jy1, Nxy, Nxy);
Jy2 = sparse(I, J, V_jy2, Nxy, Nxy);

% zero Dirichlet velocity rows (no diag=1 here)
bVel = boundaryInfo.allVelNodes(:);
if ~isempty(bVel)
    Jx1(bVel,:) = 0;
    Jx2(bVel,:) = 0;
    Jy1(bVel,:) = 0;
    Jy2(bVel,:) = 0;
end

% zeros for VP/PV blocks (convective part does not touch them)
Jx3 = spalloc(Nxy, Npr, 0);
Jy3 = spalloc(Nxy, Npr, 0);
Jc1 = spalloc(Npr, Nxy, 0);
Jc2 = spalloc(Npr, Nxy, 0);

% global convective Jacobian
Zpr = spalloc(Npr,Npr,0);
Kc  = [Jx1 Jx2 Jx3;
       Jy1 Jy2 Jy3;
       Jc1 Jc2 Zpr];
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

function pat = get_conv_pattern(elemInfo, Nxy)
% Persistent I/J for VV blocks + per-element offsets
persistent P nEl_s Nxy_s
nEl = size(elemInfo.velElements,1);
if isempty(P) || ~isequal(nEl, nEl_s) || ~isequal(Nxy, Nxy_s)
    nV  = 10; blk = nV*nV;
    I_VV = zeros(nEl*blk,1,'int32');
    J_VV = zeros(nEl*blk,1,'int32');
    off  = zeros(nEl,1,'int32');
    cursor = int32(0);
    for e=1:nEl
        Kvel = elemInfo.velElements(e,:);
        [rv,cv] = ndgrid(Kvel,Kvel);
        idx = (1:blk) + double(cursor);
        I_VV(idx) = int32(rv(:));
        J_VV(idx) = int32(cv(:));
        off(e)    = cursor;
        cursor    = cursor + int32(blk);
    end
    P = struct('I_VV',I_VV,'J_VV',J_VV,'offVV',off);
    nEl_s = nEl; Nxy_s = Nxy;
end
pat = P;
end

function y = dgemvT(A, x)
% y = A.' * x  (A: m×n, x: m×1) → y: n×1
y = (A.' * x);
end