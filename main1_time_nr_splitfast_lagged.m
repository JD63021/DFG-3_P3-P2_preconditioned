function [U, x, y, x1, y1, z1, z2, z, T, nodeInfo, elemInfo, boundaryInfo, ...
          recordedTimes, u_series, u1_series, u2_series, ...
          forceRecordedTimes, Cd_series, Cl_series] = ...
    main1_time_nr_splitfast_lagged(U1, U2, U3, D, Re, o, dt, nt, gamma, mu, rho, ...
                  nodeInfo, elemInfo, boundaryInfo, time, multiplier, force_multiplier, ...
                  corner, bcFlags, inletProfile, timeScheme, hop)
% Transient NSE with per–time-step nonlinear convergence (while-loop).
% Linear part prebuilt once; convective residual/Jacobian via fast builders:
%   - build_residual_convective_P3P2_fast
%   - build_jacobian_convective_P3P2_fast
% BDF (implicit) branch uses lagged-Jacobian Newton: reuse LU(J) across
% inner iterations; refresh only when progress stalls or every few iters.

% ==========================================================
% === Toggle mode here: 'BDF' or 'IMEX' ====================
% ==========================================================
MODE = 'BDF';   % 'IMEX' or 'BDF'

% Only override if timeScheme explicitly says IMEX
if exist('timeScheme','var') && ischar(timeScheme) && ~isempty(strtrim(timeScheme))
    ts = strtrim(timeScheme);
    if startsWith(ts,'IMEX','IgnoreCase',true)
        MODE = 'IMEX';
    end
end

% ---------------- sizes / coords ----------------
Nxy = length(nodeInfo.velocity.x);
Npr = max(elemInfo.presElements(:));

x  = nodeInfo.velocity.x;  y  = nodeInfo.velocity.y;
x1 = nodeInfo.pressure.x;  y1 = nodeInfo.pressure.y;

% ---------------- init solution & time -----------
U = [U1; U2; U3];
currentTime = time;

% history states
u_nm1 = update_bc(U, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);
u_nm2 = u_nm1;

% ---------------- Mass matrix -------------------
disp('Building mass matrix (velocity DOFs; pressure mass=0)...');
M = mass_matrix_func_s3(nodeInfo, elemInfo, boundaryInfo, rho, corner);

% ---------------- Precompute linear operator ----
% K_lin applies velocity BC rows (identity) and pressure pin
[K_lin, ~,~,~, ~,~,~, ~,~] = build_jacobian_linear_P3P2( ...
    nodeInfo, elemInfo, boundaryInfo, ...
    D, Re, o, u_nm1(1:Nxy), u_nm1(Nxy+1:2*Nxy), u_nm1(2*Nxy+1:end), ...
    gamma, mu, rho, corner);

% ---------------- recording buffers -------------
n_record        = ceil(nt / max(1,multiplier));
recordedTimes   = zeros(n_record,1);
u_series        = zeros(Nxy, n_record);
u1_series       = zeros(Nxy, n_record);
u2_series       = zeros(Nxy, n_record);
rec_idx = 1;

% ---- force sampling plan: at most ~100 samples, never every step ----
FORCE_SAMPLES     = 100;
MIN_FORCE_STRIDE  = 5;
force_stride      = max(MIN_FORCE_STRIDE, ceil(nt / FORCE_SAMPLES));
force_steps       = 1:force_stride:nt;
n_force           = numel(force_steps);
forceRecordedTimes  = zeros(n_force,1);
Cd_series           = zeros(n_force,1);
Cl_series           = zeros(n_force,1);
force_idx = 1;
prevCd = NaN; prevCl = NaN;

% ---------------- time-step control -------------
n_initial_steps = 400;
max_dt_hop      = hop;

% ---------------- Nonlinear params --------------
maxNewtonIters = 20;      % per time step
tolStep        = 1e-6;    % ||U^{k+1} - U^k||
alphaNewton    = 1.0;     % Newton relaxation
omega_hat      = 1.0;     % IMEX: relaxation for advecting field

% BDF1 warmup count (parse from timeScheme if present)
nBDF1Steps = parse_BDF1_steps(timeScheme, 1);

% ---- lagged-Jacobian controls (BDF branch) ----
lag_every     = 2;        % refresh J every k iters regardless
stall_factor  = 0.85;     % if ||ΔU|| fails to shrink by this factor, refresh J

% =======================================================
%                     TIME LOOP
% =======================================================
for it = 1:nt
    % step size
    if it <= n_initial_steps, dt1 = dt; else, dt1 = max_dt_hop; end
    currentTime = currentTime + dt1;

    % enforce BCs on old state at t^n
    u_nm1 = update_bc(u_nm1, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);

    % time discretization coefficients
    if strcmpi(MODE,'IMEX')
        if it == 1
            alpha = 1.0;                                         % BDF1
            rhs_time = (M/dt1) * u_nm1;
            Uhat = u_nm1;                                        % AB1
        else
            alpha = 3/2;                                         % BDF2
            rhs_time = (2/dt1)*M*u_nm1 - (1/(2*dt1))*M*u_nm2;
            Uhat = 2*u_nm1 - u_nm2;                              % AB2 extrapolation
        end
    else % BDF (implicit)
        if it <= nBDF1Steps
            alpha = 1.0;                                         % BDF1
            rhs_time = (M/dt1) * u_nm1;
        else
            alpha = 3/2;                                         % BDF2
            rhs_time = (2/dt1)*M*u_nm1 - (1/(2*dt1))*M*u_nm2;
        end
    end

    % ================== Nonlinear loop (per step) ==================
    U_guess = u_nm1; resNorm = Inf; iter = 0;
    prev_dU = Inf;

    % Precompute A for this step (used in both branches)
    A_step = (alpha/dt1)*M + K_lin;

    % LU cache for lagged-Jacobian (BDF branch)
    L = []; Ufac = []; p = []; q = [];  %#ok<NASGU>
    J_current = [];                      %#ok<NASGU>

    while (resNorm > tolStep && iter < maxNewtonIters)
        iter = iter + 1;

        Ug1 = U_guess(1:Nxy); Ug2 = U_guess(Nxy+1:2*Nxy); Ug3 = U_guess(2*Nxy+1:end);

        if strcmpi(MODE,'IMEX')
            % ----- IMEX Picard: explicit convection evaluated at Uhat -----
            r_conv = build_residual_convective_P3P2_fast( ...
                nodeInfo, elemInfo, boundaryInfo, ...
                D, Re, o, Uhat(1:Nxy), Uhat(Nxy+1:2*Nxy), Uhat(2*Nxy+1:end), ...
                gamma, mu, rho, currentTime, bcFlags, inletProfile, corner);

            RHS   = rhs_time - r_conv;
            U_new = A_step \ RHS;
            U_new = update_bc(U_new, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);

            resNorm = norm(U_new - U_guess);
            Uhat    = omega_hat*U_new + (1-omega_hat)*Uhat;
            U_guess = U_new;

            fprintf('  [it %4d] %s iter %2d: dt=%.3e  ||ΔU||=%.3e\n', it, MODE, iter, dt1, resNorm);

        else
            % ----- BDF Newton with lagged Jacobian -----
            % residual pieces
            r_conv = build_residual_convective_P3P2_fast( ...
                nodeInfo, elemInfo, boundaryInfo, ...
                D, Re, o, Ug1, Ug2, Ug3, gamma, mu, rho, currentTime, bcFlags, inletProfile, corner);
            R = (alpha/dt1)*(M*U_guess) + K_lin*U_guess + r_conv - rhs_time;

            % decide if we refresh K_conv & LU(J)
            need_refresh = false;
            if iter == 1
                need_refresh = true;
            elseif mod(iter, lag_every) == 0
                need_refresh = true;
            elseif isfinite(prev_dU) && resNorm > stall_factor*prev_dU
                need_refresh = true;
            end

            if need_refresh
                [K_conv, ~,~,~, ~,~,~, ~,~] = build_jacobian_convective_P3P2_fast( ...
                    nodeInfo, elemInfo, boundaryInfo, ...
                    D, Re, o, Ug1, Ug2, Ug3, gamma, mu, rho, corner);
                J = A_step + K_conv;
                [L,Ufac,p,q] = lu(J,'vector');  % sparse LU with built-in ordering
            end

            % solve J * delta = R using cached LU
            tmp   = L \ R(p);
            
            delta = zeros(size(R));
            delta(q) = Ufac \ tmp;
            

            U_new = U_guess - alphaNewton*delta;
            U_new = update_bc(U_new, boundaryInfo, nodeInfo, Nxy, currentTime, corner, bcFlags, inletProfile);

            prev_dU = resNorm;
            resNorm = norm(U_new - U_guess);
            U_guess = U_new;

            fprintf('  [it %4d] %s iter %2d: dt=%.3e  ||ΔU||=%.3e  %s\n', ...
                    it, MODE, iter, dt1, resNorm, ternary(need_refresh,'[refresh J]',''));
        end
    end

    U = U_guess;

    % shift history
    u_nm2 = u_nm1;
    u_nm1 = U;

    if any(abs(U) > 1e10)
        disp('Solution blew up; stopping.'); break;
    end

    % ---- diagnostics snapshots ----
    if mod(it, max(1,multiplier)) == 0 && rec_idx <= n_record
        recordedTimes(rec_idx) = currentTime;
        U1_final = U(1:Nxy);
        U2_final = U(Nxy+1:2*Nxy);
        Ut = hypot(U1_final, U2_final);
        u_series(:,rec_idx)  = Ut;
        u1_series(:,rec_idx) = U1_final;
        u2_series(:,rec_idx) = U2_final;
        rec_idx = rec_idx + 1;
    end

    % ---- forces: compute only at chosen steps ----
    if force_idx <= n_force && it == force_steps(force_idx)
        Uc1 = U(1:Nxy); Uc2 = U(Nxy+1:2*Nxy); Uc3 = U(2*Nxy+1:end);
        [F_x1, F_y1] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, 5, 0.001);
        [F_x2, F_y2] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, 6, 0.001);
        [F_x3, F_y3] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, 7, 0.001);
        [F_x4, F_y4] = force_int(Uc1, Uc2, Uc3, nodeInfo, elemInfo, boundaryInfo, 8, 0.001);
        F_x = F_x1 + F_x2 + F_x3 + F_x4;
        F_y = F_y1 + F_y2 + F_y3 + F_y4;

        Ua = 1; Dc = 0.1;
        Cd = F_x / (0.5 * rho * Ua^2 * Dc);
        Cl = F_y / (0.5 * rho * Ua^2 * Dc);
        prevCd = Cd; prevCl = Cl;

        forceRecordedTimes(force_idx) = currentTime;
        Cd_series(force_idx) = Cd;
        Cl_series(force_idx) = Cl;
        force_idx = force_idx + 1;

        fprintf('  ⇢ t=%.5f  ||U||=%.3e  Cd=%.3e  Cl=%.3e  [measured]\n', ...
                currentTime, norm(U), Cd, Cl);
    else
        Cd = prevCd; Cl = prevCl;
        fprintf('  ⇢ t=%.5f  ||U||=%.3e  Cd=%.3e  Cl=%.3e\n', ...
                currentTime, norm(U), Cd, Cl);
    end
end

% ---------------- Final outputs ----------------
T  = currentTime;
z1 = U(1:Nxy);
z2 = U(Nxy+1:2*Nxy);
z3 = U(2*Nxy+1:end);
z  = hypot(z1,z2);

end

% ====================== helpers =========================
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