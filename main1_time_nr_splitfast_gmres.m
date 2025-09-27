function [U, x, y, x1, y1, z1, z2, z, T, nodeInfo, elemInfo, boundaryInfo, ...
          recordedTimes, u_series, u1_series, u2_series, ...
          forceRecordedTimes, Cd_series, Cl_series] = ...
    main1_time_nr_splitfast_gmres(U1, U2, U3, D, Re, o, dt, nt, gamma, mu, rho, ...
                  nodeInfo, elemInfo, boundaryInfo, time, multiplier, force_multiplier, ...
                  corner, bcFlags, inletProfile, timeScheme, hop)
% P3–P2 NSE (triangles): Newton–GMRES with LU preconditioner that includes frozen advection.
% J = M + a*dt*(K_lin + K_conv(U)),  M1^{-1} ~ (M + a*dt*(K_lin + K_conv(Uhat)))^{-1}
% If GMRES struggles, refresh preconditioner with current U_guess and retry once.

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

% ---------------- Mass & linear operator (once) -----------
disp('Building mass matrix (velocity DOFs; pressure mass=0)...');
M = mass_matrix_func_s3(nodeInfo, elemInfo, boundaryInfo, rho, corner);

disp('Building linear operator K_lin (P3–P2) ...');
[K_lin, ~,~,~, ~,~,~, ~,~] = build_jacobian_linear_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, U_nm1(1:Nxy), U_nm1(Nxy+1:2*Nxy), U_nm1(2*Nxy+1:end), ...
    gamma, mu, rho, corner);

% ---------------- recording buffers -------------
n_record  = max(1, ceil(nt / max(1,multiplier)));
recordedTimes = zeros(n_record,1);
u_series      = zeros(Nxy, n_record);
u1_series     = zeros(Nxy, n_record);
u2_series     = zeros(Nxy, n_record);
rec_idx = 1;

FORCE_SAMPLES     = 100;
MIN_FORCE_STRIDE  = max(1, force_multiplier);
force_stride      = max(MIN_FORCE_STRIDE, ceil(nt / FORCE_SAMPLES));
force_steps       = 1:force_stride:nt;
n_force           = numel(force_steps);
forceRecordedTimes = zeros(n_force,1);
Cd_series          = zeros(n_force,1);
Cl_series          = zeros(n_force,1);
force_idx = 1;  prevCd = NaN; prevCl = NaN;

% ---------------- time-step control -------------
n_initial_steps = 40;
max_dt_hop      = hop;
nBDF1Steps      = parse_BDF1_steps(timeScheme, 1);

% ---------------- GMRES / Newton params --------
maxNewtonIters = 20;
tolStep        = 1e-6;
gmres_restart  = 50;
gmres_maxit    = 200;

% Inexact-Newton inner tol (Eisenstat–Walker-like, hard-clamped)
eta_max  = 1e-2;  eta_min = 1e-10;  ew_gamma = 0.9;  ew_beta = 1.5;

% Lag strategy for K_conv(U) (true Jacobian)
lag_every    = 1;
stall_factor = 0.85;

% Refresh preconditioner during Newton if GMRES struggles
gmres_iter_thresh_total = 60;     % compare to total iterations outer*restart + inner
gmres_relres_stall      = 1e-1;

% =======================================================
%                        TIME LOOP
% =======================================================
for it = 1:nt
    if it <= n_initial_steps, dt1 = dt; else, dt1 = max_dt_hop; end
    currentTime = currentTime + dt1;

    % --- BDF coefficients & RHS mass combo ---
    if it <= nBDF1Steps
        scheme   = 'BDF1';
        alpha    = 1.0;
        RHS_time = M * U_nm1;
        Uhat     = U_nm1;
    else
        scheme   = 'BDF2';
        alpha    = 2/3;
        RHS_time = (4/3)*M*U_nm1 - (1/3)*M*U_nm2;
        Uhat     = 2*U_nm1 - U_nm2;   % AB2 predictor
    end
    fprintf('\n=== it %d/%d | %s | dt=%.3e ===\n', it, nt, scheme, dt1);

    % --- Preconditioner: include frozen advection K_conv(Uhat) ----------
    tP = tic;
    [Kconv_hat, ~,~,~, ~,~,~, ~,~] = build_jacobian_convective_P3P2_fast( ...
        nodeInfo, elemInfo, boundaryInfo, ...
        D, Re, o, Uhat(1:Nxy), Uhat(Nxy+1:2*Nxy), Uhat(2*Nxy+1:end), ...
        gamma, mu, rho, corner);
    A_prec = M + alpha*dt1*(K_lin + Kconv_hat);
    [Lpre,Upre,pvec,qvec] = lu(A_prec,'vector');
    Mleft = @(x) apply_LU_prec(Lpre,Upre,pvec,qvec,x);
    fprintf('   [timing] prec (Kconv@Uhat + LU) in %.3fs\n', toc(tP));

    % --- Newton loop (assemble true K_conv(U) possibly lagged) ----------
    U_guess   = Uhat;       % better start than U_nm1
    stepNorm  = Inf;
    newtonIter= 0;
    prevRnorm = Inf;

    Kconv_cached = [];
    have_Kconv   = false;
    prev_dU      = Inf;

    while (stepNorm > tolStep && newtonIter < maxNewtonIters)
        newtonIter = newtonIter + 1;

        Ug1 = U_guess(1:Nxy); Ug2 = U_guess(Nxy+1:2*Nxy); Ug3 = U_guess(2*Nxy+1:end); %#ok<NASGU>

        % Residual = linear + convective
        F_lin  = K_lin * U_guess;
        F_conv = build_residual_convective_P3P2_fast( ...
                    nodeInfo, elemInfo, boundaryInfo, ...
                    D, Re, o, Ug1, Ug2, Ug3, gamma, mu, rho, currentTime, bcFlags, inletProfile, corner);
        FF     = F_lin + F_conv;

        % Dirichlet residual overwrite (inlet + walls); pin pressure
        FF = overwrite_residual_with_BCs(FF, U_guess, nodeInfo, boundaryInfo, bcFlags, inletProfile, currentTime, Nxy, corner);

        % Nonlinear residual for Newton linear solve
        R = M*U_guess + alpha*dt1*FF - RHS_time;
        Rn = norm(R);

        % Inexact-Newton target tolerance (hard clamp)
        if isfinite(prevRnorm) && prevRnorm > 0
            eta_k = ew_gamma * (Rn/prevRnorm)^ew_beta;
        else
            eta_k = eta_max;
        end
        eta_k = min(max(eta_k, eta_min), eta_max);   % clamp to [1e-10,1e-2]

        % Refresh K_conv(U) for the true Jacobian?
        need_refresh = (~have_Kconv) || (mod(newtonIter,lag_every)==0) || ...
                       (isfinite(prev_dU) && stepNorm > stall_factor*prev_dU);

        if need_refresh
            tA = tic;
            [K_conv, ~,~,~, ~,~,~, ~,~] = build_jacobian_convective_P3P2_fast( ...
                nodeInfo, elemInfo, boundaryInfo, ...
                D, Re, o, Ug1, Ug2, Ug3, gamma, mu, rho, corner);
            Kconv_cached = K_conv;
            have_Kconv   = true;
            fprintf('   [timing] K_conv(U) assembled in %.3fs\n', toc(tA));
        else
            K_conv = Kconv_cached;
        end

        % Linear operator handle J*v using assembled matrices
        Jv = @(v) ( M*v + alpha*dt1*( K_lin*v + K_conv*v ) );

        % --- GMRES solve with possible preconditioner refresh-on-fail ----
        [delta, flag, relres, iters] = gmres(Jv, R, gmres_restart, eta_k, gmres_maxit, Mleft, []);
        niter = total_iters(iters, gmres_restart);

        refreshed_once = false;
        if (flag~=0) || (niter > gmres_iter_thresh_total) || (relres > gmres_relres_stall)
            % Refresh preconditioner with current U_guess and retry once
            tPr2 = tic;
            [Kconv_hat2, ~,~,~, ~,~,~, ~,~] = build_jacobian_convective_P3P2_fast( ...
                nodeInfo, elemInfo, boundaryInfo, ...
                D, Re, o, Ug1, Ug2, Ug3, gamma, mu, rho, corner);
            A_prec2 = M + alpha*dt1*(K_lin + Kconv_hat2);
            [Lpre,Upre,pvec,qvec] = lu(A_prec2,'vector');
            Mleft = @(x) apply_LU_prec(Lpre,Upre,pvec,qvec,x);
            refreshed_once = true;
            fprintf('   [prec refresh] rebuilt at U^k in %.3fs (iters=%s → total=%d, relres=%.2e, flag=%d)\n', ...
                    toc(tPr2), vec2str(iters), niter, relres, flag);

            [delta, flag, relres, iters] = gmres(Jv, R, gmres_restart, eta_k, gmres_maxit, Mleft, []);
            niter = total_iters(iters, gmres_restart);
        end

        % Update & enforce BCs
        U_new    = U_guess - delta;
        U_new    = update_bc(U_new, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);

        prev_dU  = stepNorm;
        stepNorm = norm(U_new - U_guess);
        prevRnorm= Rn;

        fprintf('   Newton %2d: ||ΔU||=%.3e  (tol=%.1e, relres=%.2e, it=%s → total=%d)%s%s\n', ...
                newtonIter, stepNorm, eta_k, relres, vec2str(iters), niter, ...
                ternary(need_refresh,' [Kconv refresh]',''), ...
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

    % ---- forces (sample sparsely in real runs) ----
    if force_idx <= n_force && it == force_steps(force_idx)
        Uc1 = U(1:Nxy); Uc2 = U(Nxy+1:2*Nxy); Uc3 = U(2*Nxy+1:end);
        F_x = 0; F_y = 0;
        for f = [5 6 7 8]
            fn = ['flag_' num2str(f)];
            if isfield(boundaryInfo, fn) && ~isempty(boundaryInfo.(fn))
                [Fx, Fy] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, f, 1e-3);
                F_x = F_x + Fx; F_y = F_y + Fy;
            end
        end
        Ua = 1; Dc = 0.1;
        Cd = F_x/(0.5*rho*Ua*Ua*Dc);
        Cl = F_y/(0.5*rho*Ua*Ua*Dc);
        prevCd = Cd; prevCl = Cl;

        forceRecordedTimes(force_idx) = currentTime;
        Cd_series(force_idx) = Cd;
        Cl_series(force_idx) = Cl;
        force_idx = force_idx + 1;

        fprintf('  ⇢ t=%.5f  ||U||=%.3e  Cd=%.3e  Cl=%.3e  [measured]\n', ...
                currentTime, norm(U), Cd, Cl);
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

% ====================== helpers =========================
function y = apply_LU_prec(L,U,p,q,x)
% Left preconditioner solve y = (LU)^{-1} x, with lu(...,'vector')
y = zeros(size(x));
tmp = L \ x(p);
y(q) = U \ tmp;
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

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function str = vec2str(v)
if numel(v)==2, str = sprintf('[%d,%d]', v(1), v(2));
else,           str = sprintf('%d', v);
end
end

function n = total_iters(iter, restart)
% Convert GMRES iter output to a single total-iteration count.
if numel(iter) == 2
    n = iter(1)*restart + iter(2);
else
    n = iter;
end
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