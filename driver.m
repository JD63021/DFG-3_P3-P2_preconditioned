% driver (originally quadratic_solver_transient upgraded for cubic).m
% =========================================================================
% This driver sets up the mesh, initial conditions, and calls the time‐dependent 
% solver for the Navier–Stokes problem (using P3-P2 elements). It writes out a 
% separate VTK file for each recorded diagnostic snapshot and creates a PVD file 
% for time‐aware visualization.
%
% Additionally:
%   - Diagnostic snapshots (velocity, pressure, etc.) are recorded every 
%     "multiplier" time steps.
%   - Force measurements (both Cd and Cl) are recorded every "force_multiplier" 
%     time steps.
%   - Force series data (time, Cd, Cl) are appended to existing data if available,
%     ensuring a continuous record across simulation runs.
%   - Workspace snapshots are saved by the main solver at diagnostic times.
% =========================================================================

% Use all CPU threads MATLAB is allowed to use
maxNumCompThreads('automatic');           % or maxNumCompThreads(N)

% (Optional) also set BLAS/OMP env vars before heavy math starts:
setenv('OMP_NUM_THREADS', num2str(feature('numcores')));
setenv('MKL_NUM_THREADS', num2str(feature('numcores')));



tic
clc
warning off MATLAB:NonScalarInput
warning off MATLAB:nearlySingularMatrix

%% Load test case configuration
testCase = fpc();

%% User parameters
multiplier = 500;         % Diagnostic snapshot recording multiplier (multiplier 1)
force_multiplier = 1;    % Force recording multiplier (captures both Cd and Cl)
D    = 1;               % Diffusivity (or similar)
Re   = 1;
mu   = 0.001;
rho  = 1;
o    = 1;

% Check if a simulation time variable is available from a previous run.
if exist('T_sim', 'var')
    t1 = T_sim;
    fprintf('Resuming simulation from t = %.5e\n', t1);
else
    t1 = 0;
    fprintf('Starting simulation from t = 0\n');
end

%% Generate Mesh & Setup
[nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh_p3(...
    testCase.gmshVelFile, testCase.gmshPresFile, ...
    testCase.excludedVelFlags, testCase.pBoundFlags, testCase.pBoundVals);
Nxy = length(nodeInfo.velocity.x);  % Number of velocity DOFs
Npr = max(elemInfo.presElements(:));
% For Q2–Q1 elements, the "corner" parameter may be unused:
% corner = testCase.corner;  % (if needed)

% Set initial conditions.
if exist('U','var')
    U1 = U(1:Nxy); 
    U2 = U(Nxy+1:2*Nxy); 
    U3 = U(2*Nxy+1:end);
    fprintf('Resuming simulation with existing solution U\n');
else
    U1 = zeros(Nxy,1);
    U2 = zeros(Nxy,1);
    U3 = zeros(Npr,1);
    if isfield(boundaryInfo, ['flag_' num2str(testCase.boundaryFlags.inlet(1))])
        bV = boundaryInfo.(['flag_' num2str(testCase.boundaryFlags.inlet(1))]);
        U1(bV) = 0;
    end
    fprintf('Starting simulation with initial conditions\n');
end

%% Time-stepping parameters
dt    = 0.5e-2;
hop   = 0.5e-2;
nt    = 1600;
gamma = 0;

%% Call the time-stepping solver.
% [U, x, y, x1, y1, z1, z2, z, T, nodeInfo, elemInfo, boundaryInfo, ...
%           recordedTimes, u_series, u1_series, u2_series, ...
%           Cl_recordedTimes, Cl_series] = ...
%     main1_time_nr_splitfast_gmres(U1, U2, U3, D, Re, o, dt, nt, gamma, mu, rho, ...
%                   nodeInfo, elemInfo, boundaryInfo, t1, multiplier, force_multiplier, ...
%                   testCase.corner, testCase.boundaryFlags, testCase.inletProfile, 'BDF', hop);
[U, x, y, x1, y1, z1, z2, z, T, nodeInfo, elemInfo, boundaryInfo, ...
  recordedTimes, u_series, u1_series, u2_series, ...
  forceRecordedTimes, Cd_series, Cl_series] = ...
  main1_time_nr_splitfast_gmres(U1, U2, U3, D, Re, o, dt, nt, gamma, mu, rho, ...
    nodeInfo, elemInfo, boundaryInfo, t1, multiplier, force_multiplier, ...
    testCase.corner, testCase.boundaryFlags, testCase.inletProfile, 'BDF', hop);

              
% Update velocity and pressure parts from the solution.
U1 = U(1:Nxy);
U2 = U(Nxy+1:2*Nxy);
U3 = U(2*Nxy+1:end);
u = sqrt(U1.^2+ U2.^2);
T_sim = T;
fprintf('Simulation ended at time T = %.5e\n', T_sim);

%% --- Write VTP Files for Diagnostic Snapshots ---
numSnapshots = length(recordedTimes);
for i = 1:numSnapshots
    filename = sprintf('solution_%04d.vtp', i);
    writeVTP(x, y, u_series(:, i), filename);
    fprintf('Wrote snapshot %d to file %s (time = %f)\n', i, filename, recordedTimes(i));
end

%% --- Generate a PVD File Linking the Snapshots ---
pvdFilename = 'solutions.pvd';
baseFilename = 'solution_%04d.vtp';
writePVD(recordedTimes, baseFilename, pvdFilename);

%% --- Plot the Final Velocity and Pressure Fields ---
figure;
scatter3(x, y, u, 20, u, 'filled');
title('Final Velocity Magnitude');
grid on; colormap(jet); colorbar;

figure;
scatter3(x1, y1, U3, 20, U3, 'filled');
title('Final Pressure Field');
grid on; colormap(jet); colorbar;

%% --- Append and Save Force/Lift Series ---
% Keep only actually recorded rows (nonzero time & finite values)
mask = (forceRecordedTimes > 0) & isfinite(Cd_series) & isfinite(Cl_series);
tF = forceRecordedTimes(mask);
CD = Cd_series(mask);
CL = Cl_series(mask);

% Nothing recorded? create empty 0x3 to avoid concat errors
if isempty(tF), force_data = zeros(0,3);
else,           force_data = [tF(:), CD(:), CL(:)];  % [time, Cd, Cl]
end

fname = 'force_series.mat';
if exist(fname, 'file')
    S = load(fname);
    if isfield(S, 'force_total_series') && ~isempty(S.force_total_series)
        F = S.force_total_series;
        % Ensure 3 columns: [time, Cd, Cl]
        if size(F,2) < 3, F = [F, nan(size(F,1), 3 - size(F,2))]; end
        if size(force_data,2) < 3
            force_data = [force_data, nan(size(force_data,1), 3 - size(force_data,2))];
        end
        force_total_series = [F; force_data];
    else
        force_total_series = force_data;
    end
else
    force_total_series = force_data;
end

save(fname, 'force_total_series');
fprintf('Force series saved to %s (%d rows, %d cols).\n', fname, size(force_total_series,1), size(force_total_series,2));

toc

