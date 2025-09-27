function [U, x, y, x1, y1, z1, z2, z, T, nodeInfo, elemInfo, boundaryInfo, ...
          recordedTimes, u_series, u1_series, u2_series, ...
          forceRecordedTimes, Cd_series, Cl_series] = ...
    main1_time_nr_p3p2_matrixfree(U1, U2, U3, D, Re, o, dt, nt, gamma, mu, rho, ...
                  nodeInfo, elemInfo, boundaryInfo, time, multiplier, force_multiplier, ...
                  corner, bcFlags, inletProfile, timeScheme, hop)
% P3–P2 NSE on triangles: Newton–GMRES with matrix-free convective Jacobian.
% Preconditioner includes a FROZEN ADVECTIVE matrix K_adv(Uhat); we reuse its LU
% across steps if the predictor changes little, and refresh if GMRES struggles.

% ---------------- sizes / coords ----------------
Nxy = length(nodeInfo.velocity.x);
Npr = max(elemInfo.presElements(:));
NN  = 2*Nxy + Npr;

x  = nodeInfo.velocity.x;  y  = nodeInfo.velocity.y;
x1 = nodeInfo.pressure.x;  y1 = nodeInfo.pressure.y;

% ---------------- init solution & time -----------
U = [U1; U2; U3];
currentTime = time;

% history (enforce BCs at start)
U_nm1 = update_bc(U, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);
U_nm2 = U_nm1;

% ---------------- Mass & linear operator (once) --------
disp('Building mass matrix (velocity DOFs; pressure mass=0)...');
M = mass_matrix_func_s3(nodeInfo, elemInfo, boundaryInfo, rho, corner);

disp('Building linear operator K_lin (P3–P2) ...');
[K_lin, ~,~,~, ~,~,~, ~,~] = build_jacobian_linear_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U_nm1(1:Nxy), U_nm1(Nxy+1:2*Nxy), U_nm1(2*Nxy+1:end), ...
    gamma, mu, rho, corner);

% ---------------- Convective geometry cache -------------
disp('Caching convective geometry (P3–P2 triangles) ...');
convCache = p3p2_make_convective_cache(nodeInfo, elemInfo);

% ---------------- recording buffers ---------------------
n_record  = max(1, ceil(nt / max(1,multiplier)));
recordedTimes = zeros(n_record,1);
u_series      = zeros(Nxy, n_record);
u1_series     = zeros(Nxy, n_record);
u2_series     = zeros(Nxy, n_record);
rec_idx = 1;

% FORCE_SAMPLES     = 500;
% MIN_FORCE_STRIDE  = max(1, force_multiplier);
% --- force sampling schedule (single knob = force_multiplier from driver) ---
force_stride      = max(1, force_multiplier);   % sample every 'force_multiplier' steps
force_steps       = 1:force_stride:nt;
n_force           = numel(force_steps);
forceRecordedTimes = nan(n_force,1);
Cd_series          = nan(n_force,1);
Cl_series          = nan(n_force,1);
force_idx = 1;  prevCd = NaN; prevCl = NaN;


% ---------------- time-step control ---------------------
n_initial_steps = 40;
max_dt_hop      = hop;
nBDF1Steps      = parse_BDF1_steps(timeScheme, 1);

% ---------------- GMRES / Newton params ----------------
maxNewtonIters = 20;
tolStep        = 1e-6;          % ||U^{k+1} - U^k||
gmres_restart  = 50;
gmres_maxit    = 200;

% Inexact-Newton target for GMRES (clamped)
eta_max  = 1e-2;  eta_min = 1e-10;  ew_gamma = 0.9;  ew_beta = 1.5;

% Preconditioner reuse / rebuild policy
Uhat_reuse_tol = 1e-3;          % relative change threshold
gmres_iter_thresh_total = 120;  % if more than this, refresh prec
gmres_relres_stall      = 1e-1;

% Preconditioner state (persist across steps)
Lpre=[]; Upre=[]; pvec=[]; qvec=[];
Kadv_hat = []; last_Uhat = [];

% =======================================================
%                        TIME LOOP
% =======================================================
for it = 1:nt
    % step size
    if it <= n_initial_steps, dt1 = dt; else, dt1 = max_dt_hop; end
    currentTime = currentTime + dt1;

    % scheme & RHS mass combo + AB predictor
    if it <= nBDF1Steps
        scheme   = 'BDF1';
        alpha    = 1.0;
        RHS_time = M * U_nm1;
        Uhat     = U_nm1;
    else
        scheme   = 'BDF2';
        alpha    = 2/3;
        RHS_time = (4/3)*M*U_nm1 - (1/3)*M*U_nm2;
        Uhat     = 2*U_nm1 - U_nm2;  % AB2 predictor
    end
    fprintf('\n=== it %d/%d | %s | dt=%.3e ===\n', it, nt, scheme, dt1);

    % --------- Preconditioner: K_adv(Uhat) + reuse across steps ----------
    need_prec = isempty(Lpre) || isempty(last_Uhat);
    if ~need_prec
        dUhat = norm(Uhat - last_Uhat);
        need_prec = (dUhat > Uhat_reuse_tol * max(1e-14, norm(last_Uhat)));
    end
    if need_prec
        tP = tic;
        Kadv_hat = build_advective_matrix_p3p2_cached(convCache, Re, rho, ...
                        Uhat(1:Nxy), Uhat(Nxy+1:2*Nxy), boundaryInfo);
        Kadv_full = blkdiag(Kadv_hat, sparse(Npr,Npr));
        A_prec = M + alpha*dt1*(K_lin + Kadv_full);
        [Lpre,Upre,pvec,qvec] = lu(A_prec,'vector');
        last_Uhat = Uhat;
        fprintf('   [timing] prec (Kadv@Uhat + LU) in %.3fs\n', toc(tP));
    else
        fprintf('   [reuse] preconditioner from previous step\n');
    end
    Mleft = @(x) apply_LU_prec(Lpre,Upre,pvec,qvec,x);

    % --------- Newton loop (matrix-free convective J) ----------
    U_guess   = Uhat;      % good start
    stepNorm  = Inf;
    newtonIter= 0;
    prevRnorm = Inf;

    while (stepNorm > tolStep && newtonIter < maxNewtonIters)
        newtonIter = newtonIter + 1;

        % Residual parts
        F_lin  = K_lin * U_guess;
        F_conv = build_residual_convective_p3p2_cached(convCache, Re, rho, ...
                        U_guess(1:Nxy), U_guess(Nxy+1:2*Nxy), Npr);
        FF = F_lin + F_conv;

        % Dirichlet overwrite + pinned pressure
        FF = overwrite_residual_with_BCs(FF, U_guess, nodeInfo, boundaryInfo, bcFlags, inletProfile, currentTime, Nxy, corner);

        % Nonlinear residual for Newton linear solve
        R = M*U_guess + alpha*dt1*FF - RHS_time;
        Rn = norm(R);

        % Inexact-Newton tolerance
        if isfinite(prevRnorm) && prevRnorm > 0
            eta_k = ew_gamma * (Rn/prevRnorm)^ew_beta;
        else
            eta_k = eta_max;
        end
        eta_k = min(max(eta_k, eta_min), eta_max);

        % Matrix-free J*v
        data = struct('M',M,'Klin',K_lin,'alpha',alpha,'dt1',dt1, ...
                      'cache',convCache,'Re',Re,'rho',rho,'U',U_guess, ...
                      'boundaryInfo',boundaryInfo,'corner',corner,'Nxy',Nxy,'NN',NN);
        Jv = @(v) ( data.M*v + data.alpha*data.dt1*( data.Klin*v + ...
               jconv_times_v_p3p2_cached(data.cache, data.Re, data.rho, data.U, v, data.boundaryInfo, data.corner) ) );

        % GMRES solve
        [delta, flag, relres, iters] = gmres(Jv, R, gmres_restart, eta_k, gmres_maxit, Mleft, []);
        niter = total_iters(iters, gmres_restart);

        % If GMRES struggled, refresh preconditioner at U_guess and retry once
        refreshed_once = false;
        if (flag~=0) || (niter > gmres_iter_thresh_total) || (relres > gmres_relres_stall)
            tPr2 = tic;
            Kadv_now = build_advective_matrix_p3p2_cached(convCache, Re, rho, ...
                            U_guess(1:Nxy), U_guess(Nxy+1:2*Nxy), boundaryInfo);
            Kadv_full2 = blkdiag(Kadv_now, sparse(Npr,Npr));
            A_prec2 = M + alpha*dt1*(K_lin + Kadv_full2);
            [Lpre,Upre,pvec,qvec] = lu(A_prec2,'vector');
            Mleft = @(x) apply_LU_prec(Lpre,Upre,pvec,qvec,x);
            refreshed_once = true;
            last_Uhat = U_guess;   % advance reuse anchor
            fprintf('   [prec refresh] rebuilt at U^k in %.3fs (iters=%s → total=%d, relres=%.2e, flag=%d)\n', ...
                    toc(tPr2), vec2str(iters), niter, relres, flag);

            [delta, flag, relres, iters] = gmres(Jv, R, gmres_restart, eta_k, gmres_maxit, Mleft, []);
            niter = total_iters(iters, gmres_restart);
        end

        % Update & enforce BCs
        U_new    = U_guess - delta;
        U_new    = update_bc(U_new, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);

        stepNorm  = norm(U_new - U_guess);
        prevRnorm = Rn;

        fprintf('   Newton %2d: ||ΔU||=%.3e  (tol=%.1e, relres=%.2e, it=%s → total=%d)%s\n', ...
                newtonIter, stepNorm, eta_k, relres, vec2str(iters), niter, ...
                ternary(refreshed_once,' [prec refresh]',''));

        U_guess = U_new;
    end

    % accept
    U = U_guess;

    % shift history
    U_nm2 = U_nm1;
    U_nm1 = U;

    if any(~isfinite(U)) || any(abs(U) > 1e10)
        warning('Solution diverged, stopping.'); break;
    end

    % ---- diagnostics snapshots ----
    if mod(it, max(1,multiplier)) == 0 && rec_idx <= n_record
        recordedTimes(rec_idx) = currentTime;
        U1f = U(1:Nxy); U2f = U(Nxy+1:2*Nxy);
        Ut  = hypot(U1f, U2f);
        u_series(:,rec_idx)  = Ut;
        u1_series(:,rec_idx) = U1f;
        u2_series(:,rec_idx) = U2f;
        rec_idx = rec_idx + 1;
    end

    % ---- forces (optional; keep lightweight sampling) ----
   if force_idx <= n_force && it == force_steps(force_idx)
    Uc1 = U(1:Nxy); Uc2 = U(Nxy+1:2*Nxy); Uc3 = U(2*Nxy+1:end);

    % Sum over cylinder flags that exist
    F_x = 0; F_y = 0;
    for f = [5 6 7 8]
        fn = ['flag_' num2str(f)];
        haveLines = isfield(boundaryInfo,'velLine4Elements') && ...
                    isfield(boundaryInfo.velLine4Elements, fn) && ...
                    ~isempty(boundaryInfo.velLine4Elements.(fn));
        if haveLines && exist('force_int','file')==2
            % IMPORTANT: pass the same cRef you used before (see C below)
            [Fx, Fy] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, f, mu, [0.2, 0.2]);
            F_x = F_x + Fx; F_y = F_y + Fy;
        end
    end

    Ua = 1; Dc = 0.1;
    Cd = F_x/(0.5*rho*Ua*Ua*Dc);
    Cl = F_y/(0.5*rho*Ua*Ua*Dc);

    % STORE FIRST
    forceRecordedTimes(force_idx) = currentTime;
    Cd_series(force_idx) = Cd;
    Cl_series(force_idx) = Cl;
    prevCd = Cd; prevCl = Cl;
    force_idx = force_idx + 1;

    % PRINT *FROM THE STORED VALUES*
    fprintf('  ⇢ t=%.5f  ||U||=%.3e  Cd=%.3e  Cl=%.3e  [measured]\n', ...
            currentTime, norm(U), Cd_series(force_idx-1), Cl_series(force_idx-1));
else
    fprintf('  ⇢ t=%.5f  ||U||=%.3e  Cd=%s  Cl=%s\n', ...
            currentTime, norm(U), num2str(prevCd,'%.3e'), num2str(prevCl,'%.3e'));
end


end

% ---------------- Final outputs ----------------
T  = currentTime;
z1 = U(1:Nxy);
z2 = U(Nxy+1:2*Nxy);
z  = hypot(z1,z2);
end % ====== main ======

% ======================================================================
%                             HELPERS
% ======================================================================

function y = jconv_times_v_p3p2_cached(cache, Re, rho, U, v, boundaryInfo, corner)
% Matrix-free action of the convective Jacobian on v (P3–P2 triangles)
Nxy = cache.Nxy; Npr = cache.Npr;
y_u = zeros(Nxy,1);
y_v = zeros(Nxy,1);

U1 = U(1:Nxy); U2 = U(Nxy+1:2*Nxy);
v1 = v(1:Nxy); v2 = v(Nxy+1:2*Nxy);

for e = 1:cache.numEl
    K   = cache.el(e).Kvel;
    Ni  = cache.el(e).Ni;     % 10 x G
    dNx = cache.el(e).dNx;    % 10 x G
    dNy = cache.el(e).dNy;    % 10 x G
    W   = cache.el(e).W(:);   % G x 1

    U1e = U1(K); U2e = U2(K);
    v1e = v1(K); v2e = v2(K);

    a1   = Ni.'  * U1e;  a2   = Ni.'  * U2e;   % Gx1
    du1x = dNx.' * U1e;  du1y = dNy.' * U1e;
    du2x = dNx.' * U2e;  du2y = dNy.' * U2e;

    b1   = Ni.'  * v1e;  b2   = Ni.'  * v2e;
    dv1x = dNx.' * v1e;  dv1y = dNy.' * v1e;
    dv2x = dNx.' * v2e;  dv2y = dNy.' * v2e;

    cx = a1.*dv1x + a2.*dv1y + b1.*du1x + b2.*du1y;
    cy = a1.*dv2x + a2.*dv2y + b1.*du2x + b2.*du2y;

    y_u(K) = y_u(K) + rho*Re * (Ni * (cx .* W));
    y_v(K) = y_v(K) + rho*Re * (Ni * (cy .* W));
end

y = [y_u; y_v; zeros(Npr,1)];

% Dirichlet rows & pinned pressure row → identity in J
if isfield(boundaryInfo,'allVelNodes') && ~isempty(boundaryInfo.allVelNodes)
    b = boundaryInfo.allVelNodes(:);
    b = b(b>=1 & b<=Nxy);
    y(b)         = v(b);
    y(Nxy + b)   = v(Nxy + b);
end
y(2*cache.Nxy + corner) = v(2*cache.Nxy + corner);
end

function F_conv = build_residual_convective_p3p2_cached(cache, Re, rho, U1, U2, Npr)
% Convective residual ONLY, using cached geometry (P3–P2)
Nxy = cache.Nxy;
V = zeros(Nxy,1);
W = zeros(Nxy,1);
for e = 1:cache.numEl
    K   = cache.el(e).Kvel;
    Ni  = cache.el(e).Ni;     % 10 x G
    dNx = cache.el(e).dNx;    % 10 x G
    dNy = cache.el(e).dNy;    % 10 x G
    WG  = cache.el(e).W(:);   % G x 1

    U1e = U1(K); U2e = U2(K);

    a1   = Ni.'  * U1e;  a2   = Ni.'  * U2e;   % Gx1
    du1x = dNx.' * U1e;  du1y = dNy.' * U1e;
    du2x = dNx.' * U2e;  du2y = dNy.' * U2e;

    cu = a1.*du1x + a2.*du1y;
    cv = a1.*du2x + a2.*du2y;

    V(K) = V(K) + rho*Re * (Ni * (cu .* WG));
    W(K) = W(K) + rho*Re * (Ni * (cv .* WG));
end
F_conv = [V; W; zeros(Npr,1)];
end

function Kadv = build_advective_matrix_p3p2_cached(cache, Re, rho, U1, U2, boundaryInfo)
% Assemble velocity-block (2Nxy x 2Nxy) matrix for (u·∇)(·), using cached Ni,dNx,dNy,W
Nxy = cache.Nxy; nEl = cache.numEl;
nzPerEl = 10*10;

I = zeros(2*nzPerEl*nEl,1);
J = zeros(2*nzPerEl*nEl,1);
S = zeros(2*nzPerEl*nEl,1);
pos = 0;

for e = 1:nEl
    K   = cache.el(e).Kvel;
    Ni  = cache.el(e).Ni;   % 10 x G
    dNx = cache.el(e).dNx;  % 10 x G
    dNy = cache.el(e).dNy;  % 10 x G
    WG  = cache.el(e).W(:); % G x 1

    U1e = U1(K); U2e = U2(K);
    a1  = Ni.' * U1e;       % G x 1
    a2  = Ni.' * U2e;

    Gadv = zeros(10,10);
    for g = 1:numel(WG)
        gradT = (a1(g)*dNx(:,g) + a2(g)*dNy(:,g)).';  % 1 x 10
        Gadv  = Gadv + rho*Re * WG(g) * (Ni(:,g) * gradT);  % 10x10
    end

    [rows, cols] = ndgrid(K, K);
    nn = numel(Gadv);
    I(pos+(1:nn)) = rows(:);
    J(pos+(1:nn)) = cols(:);
    S(pos+(1:nn)) = Gadv(:);
    pos = pos + nn;

    rows = rows + Nxy; cols = cols + Nxy;
    I(pos+(1:nn)) = rows(:);
    J(pos+(1:nn)) = cols(:);
    S(pos+(1:nn)) = Gadv(:);
    pos = pos + nn;
end

Kadv = sparse(I(1:pos), J(1:pos), S(1:pos), 2*Nxy, 2*Nxy);

% Zero Dirichlet velocity rows to preserve identity from K_lin in A_prec
if isfield(boundaryInfo,'allVelNodes') && ~isempty(boundaryInfo.allVelNodes)
    b = boundaryInfo.allVelNodes(:);
    b = b(b>=1 & b<=Nxy);
    Kadv(b,:)       = 0;
    Kadv(Nxy + b,:) = 0;
end
end

function cache = p3p2_make_convective_cache(nodeInfo, elemInfo)
% Cache per-element Ni, dNx, dNy and W for P3 triangles with Dunavant-7 rule
Nxy = length(nodeInfo.velocity.x);
Npr = max(elemInfo.presElements(:));
nEl = size(elemInfo.velElements,1);

[lam, w] = dunavant7();   % G x 3 barycentric and G x 1 weights
G = size(lam,1);

cache.el(nEl) = struct('Ni',[],'dNx',[],'dNy',[],'W',[],'Kvel',[]);
for e = 1:nEl
    Kvel = elemInfo.velElements(e,:);   % 10 P3 nodes
    xE   = nodeInfo.velocity.x(Kvel);
    yE   = nodeInfo.velocity.y(Kvel);

    Ni  = zeros(10, G);
    dNx = zeros(10, G);
    dNy = zeros(10, G);
    Wg  = zeros(G, 1);

    for g = 1:G
        L1 = lam(g,1); L2 = lam(g,2); L3 = lam(g,3);
        xi = L2; eta = L3;  % reference tri: (0,0),(1,0),(0,1)

        [Nref, dNxi, dNeta] = p3_basis_ref(xi, eta);

        % Jacobian at gp
        dX_dxi  = sum(xE .* dNxi);
        dX_deta = sum(xE .* dNeta);
        dY_dxi  = sum(yE .* dNxi);
        dY_deta = sum(yE .* dNeta);

        detJ = dX_dxi*dY_deta - dX_deta*dY_dxi;
        invJ = [ dY_deta, -dY_dxi; -dX_deta, dX_dxi ] / detJ;

        dNx(:,g) = invJ(1,1)*dNxi + invJ(1,2)*dNeta;
        dNy(:,g) = invJ(2,1)*dNxi + invJ(2,2)*dNeta;

        Ni(:,g)  = Nref;
        Wg(g)    = detJ * w(g);
    end

    cache.el(e).Ni   = Ni;
    cache.el(e).dNx  = dNx;
    cache.el(e).dNy  = dNy;
    cache.el(e).W    = Wg;
    cache.el(e).Kvel = Kvel(:);
end
cache.Nxy   = Nxy;
cache.Npr   = Npr;
cache.numEl = nEl;
end

function [lam, w] = dunavant7()
% 7-point Dunavant rule on reference triangle (area 1/2)
% Returns barycentric coordinates lam(:,1:3) and weights w(:)
lam = [ 1/3, 1/3, 1/3;
        0.0597158717898, 0.4701420641051, 0.4701420641051;
        0.4701420641051, 0.0597158717898, 0.4701420641051;
        0.4701420641051, 0.4701420641051, 0.0597158717898;
        0.7974269853531, 0.1012865073235, 0.1012865073235;
        0.1012865073235, 0.7974269853531, 0.1012865073235;
        0.1012865073235, 0.1012865073235, 0.7974269853531 ];
w = [ 0.2250000000000;
      0.1323941527885;
      0.1323941527885;
      0.1323941527885;
      0.1259391805448;
      0.1259391805448;
      0.1259391805448 ];
end

function [N, dNxi, dNeta] = p3_basis_ref(xi, eta)
% P3 triangle basis on reference tri with vertices (0,0),(1,0),(0,1)
L1 = 1 - xi - eta;
L2 = xi;
L3 = eta;

N = zeros(10,1);
% corners
N(1) = L1*(3*L1 - 1)*(3*L1 - 2)/2;
N(2) = L2*(3*L2 - 1)*(3*L2 - 2)/2;
N(3) = L3*(3*L3 - 1)*(3*L3 - 2)/2;
% mid-edges (two per edge) + bubble
c = 9/2;
N(4) = c * L1*L2*(3*L1 - 1);
N(5) = c * L1*L2*(3*L2 - 1);
N(6) = c * L2*L3*(3*L2 - 1);
N(7) = c * L2*L3*(3*L3 - 1);
N(8) = c * L3*L1*(3*L3 - 1);
N(9) = c * L3*L1*(3*L1 - 1);
N(10)= 27*L1*L2*L3;

dL1_dxi=-1; dL1_deta=-1;
dL2_dxi= 1; dL2_deta= 0;
dL3_dxi= 0; dL3_deta= 1;

dNxi  = zeros(10,1); dNeta = zeros(10,1);
% corner1
f=(3*L1 -1)*(3*L1 -2)/2; df_dL1=(9*(2*L1 -1))/2;
dNxi(1)  = dL1_dxi*f + L1*(df_dL1*dL1_dxi);
dNeta(1) = dL1_deta*f + L1*(df_dL1*dL1_deta);
% corner2
f2=(3*L2 -1)*(3*L2 -2)/2; df2=(9*(2*L2 -1))/2;
dNxi(2)  = dL2_dxi*f2 + L2*(df2*dL2_dxi);
dNeta(2) = dL2_deta*f2 + L2*(df2*dL2_deta);
% corner3
f3=(3*L3 -1)*(3*L3 -2)/2; df3=(9*(2*L3 -1))/2;
dNxi(3)  = dL3_dxi*f3 + L3*(df3*dL3_dxi);
dNeta(3) = dL3_deta*f3 + L3*(df3*dL3_deta);
% mid/bubble via product rule helpers
[dNxi(4), dNeta(4)] = d_of_prod3(c, L1, L2, (3*L1 -1), ...
                                 dL1_dxi, dL1_deta, dL2_dxi, dL2_deta, 3*dL1_dxi, 3*dL1_deta);
[dNxi(5), dNeta(5)] = d_of_prod3(c, L1, L2, (3*L2 -1), ...
                                 dL1_dxi, dL1_deta, dL2_dxi, dL2_deta, 3*dL2_dxi, 3*dL2_deta);
[dNxi(6), dNeta(6)] = d_of_prod3(c, L2, L3, (3*L2 -1), ...
                                 dL2_dxi, dL2_deta, dL3_dxi, dL3_deta, 3*dL2_dxi, 3*dL2_deta);
[dNxi(7), dNeta(7)] = d_of_prod3(c, L2, L3, (3*L3 -1), ...
                                 dL2_dxi, dL2_deta, dL3_dxi, dL3_deta, 3*dL3_dxi, 3*dL3_deta);
[dNxi(8), dNeta(8)] = d_of_prod3(c, L3, L1, (3*L3 -1), ...
                                 dL3_dxi, dL3_deta, dL1_dxi, dL1_deta, 3*dL3_dxi, 3*dL3_deta);
[dNxi(9), dNeta(9)] = d_of_prod3(c, L3, L1, (3*L1 -1), ...
                                 dL3_dxi, dL3_deta, dL1_dxi, dL1_deta, 3*dL1_dxi, 3*dL1_deta);
[dNxi(10),dNeta(10)] = d_of_prod3(27, L1, L2, L3, ...
                                  dL1_dxi, dL1_deta, dL2_dxi, dL2_deta, dL3_dxi, dL3_deta);
end

function [dNdxi, dNeta] = d_of_prod3(scale, X, Y, Z, dXdxi, dXdeta, dYdxi, dYdeta, dZdxi, dZdeta)
dNdxi = scale*( dXdxi*Y*Z + X*dYdxi*Z + X*Y*dZdxi );
dNeta = scale*( dXdeta*Y*Z + X*dYdeta*Z + X*Y*dZdeta );
end

function y = apply_LU_prec(L,U,p,q,x)
% Left preconditioner solve y = (LU)^{-1} x, with lu(...,'vector')
y = zeros(size(x));
tmp = L \ x(p);
y(q) = U \ tmp;
end

function FF = overwrite_residual_with_BCs(FF, U_guess, nodeInfo, boundaryInfo, bcFlags, inletProfile, t, Nxy, corner)
% Dirichlet velocity + pinned pressure row in residual
% inlet
if isfield(bcFlags,'inlet')
    topFlag = bcFlags.inlet;
    if isfield(boundaryInfo, ['flag_' num2str(topFlag)])
        inletNodes = boundaryInfo.(['flag_' num2str(topFlag)])(:);
    else
        inletNodes = [];
    end
else
    inletNodes = [];
end
if ~isempty(inletNodes)
    yMin = min(nodeInfo.velocity.y);
    yMax = max(nodeInfo.velocity.y);
    H    = yMax - yMin;
    yvals = nodeInfo.velocity.y(inletNodes);
    Uin  = arrayfun(@(yy) inletProfile(t, yy, H), yvals);
    U1blk = U_guess(1:Nxy);
    U2blk = U_guess(Nxy+1:2*Nxy);
    FF(inletNodes)       = U1blk(inletNodes) - Uin(:);
    FF(Nxy + inletNodes) = U2blk(inletNodes);
end
% walls / sides
sideNodes = [];
if isfield(bcFlags,'wall')
    for k = 1:length(bcFlags.wall)
        fn = ['flag_' num2str(bcFlags.wall(k))];
        if isfield(boundaryInfo, fn)
            sideNodes = [sideNodes; boundaryInfo.(fn)(:)]; %#ok<AGROW>
        end
    end
    sideNodes = unique(sideNodes);
end
if ~isempty(sideNodes)
    U1blk = U_guess(1:Nxy);
    U2blk = U_guess(Nxy+1:2*Nxy);
    FF(sideNodes)       = U1blk(sideNodes);
    FF(Nxy + sideNodes) = U2blk(sideNodes);
end
% pinned pressure
FF(2*Nxy + corner) = U_guess(2*Nxy + corner);
end

function nBDF1Steps = parse_BDF1_steps(timeScheme, defaultVal)
nBDF1Steps = defaultVal;
if isempty(timeScheme) || ~ischar(timeScheme), return; end
tok = regexp(timeScheme,'(?i)bdf1\s*[\(\:x]\s*(\d+)\s*[\)]?','tokens','once');
if ~isempty(tok), nBDF1Steps = max(1, str2double(tok{1})); end
if strcmpi(strtrim(timeScheme),'BDF1')
    nBDF1Steps = intmax; % never switch to BDF2
end
end

function str = vec2str(v)
if numel(v)==2, str = sprintf('[%d,%d]', v(1), v(2));
else,           str = sprintf('%d', v);
end
end

function n = total_iters(iter, restart)
if numel(iter) == 2
    n = iter(1)*restart + iter(2);
else
    n = iter;
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
