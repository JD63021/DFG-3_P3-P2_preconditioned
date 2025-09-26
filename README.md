# P3–P2 Navier–Stokes (2D, Triangles) — Matrix-Free Newton/GMRES

High-order incompressible Navier–Stokes solver on **P3 velocity / P2 pressure** (10-node / 6-node triangles).
Time integration uses BDF1/BDF2 with **matrix-free** convective Jacobian in GMRES and an **LU preconditioner** that includes a *frozen* advective operator.

The code writes Paraview-ready VTP snapshots and a PVD time index, and (optionally) drag/lift time series.

---

## Highlights

* **Elements:** P3 (velocity) / P2 (pressure) on triangles.
* **Newton/GMRES:** Inexact Newton (Eisenstat–Walker) + restarted GMRES.
* **Preconditioner:** $M + \alpha\Delta t\,(K_{\text{lin}} + K_{\text{conv}}( \hat u ))$ with sparse LU.
* **Convective terms:** matrix-free J·v, fast assembled residual/Jacobian for refreshes.
* **BCs:** inlet profile, walls/no-slip, pinned pressure DOF.
* **I/O:** `solution_####.vtp` + `solutions.pvd`, plus optional `force_series.mat`.

---

## Requirements

* MATLAB R2020b+ (tested on recent releases).
* Gmsh meshes exported to MATLAB `.m` structs:

  * **Velocity mesh:** `TRIANGLES10` + `LINES4`
  * **Pressure mesh:** `TRIANGLES6` + `LINES3`

No external toolboxes are required.

---

## Repo Layout (key files)

```
driver.m                         % top-level script (run this)
main1_time_nr_splitfast_gmres.m  % P3–P2 time integrator (Newton/GMRES)
mass_matrix_func_s3.m            % velocity mass (pressure mass = 0)
build_jacobian_linear_P3P2.m     % viscous + pressure/continuity (P3–P2)
build_residual_convective_P3P2_fast.m
build_jacobian_convective_P3P2_fast.m
mesh5_gmsh_p3.m                  % Gmsh reader → nodeInfo/elemInfo/boundaryInfo
writeVTP.m, writePVD.m           % Paraview output
force_int.m                      % (optional) boundary traction integration
fpc.m                            % case configuration (file paths, flags, inlet)
```

> If you keep a different file organization, adjust `driver.m` paths accordingly.

---

## Quick Start

1. **Point to your meshes** in `fpc.m`:

   ```matlab
   function cfg = fpc()
     cfg.gmshVelFile   = 'mesh_vel_p3.m';   % has TRIANGLES10, LINES4
     cfg.gmshPresFile  = 'mesh_pre_p2.m';   % has TRIANGLES6, LINES3
     cfg.excludedVelFlags = [];             % optional
     cfg.pBoundFlags   = []; cfg.pBoundVals = [];  % optional pressure BCs

     % Boundary flag ids from Gmsh
     cfg.boundaryFlags.inlet = 1;
     cfg.boundaryFlags.wall  = [2 3 4];

     % Inlet profile: function handle @(t,y,H) → u_in(y)
     cfg.inletProfile = @(t,y,H) 4*(y.*(H - y))/H^2;  % parabola example

     % Pinned pressure DOF index (in P2 space)
     cfg.corner = 1;
   end
   ```

2. **Run**:

   ```matlab
   >> driver
   ```

   You’ll see per-step logs, GMRES stats, and files `solution_####.vtp` + `solutions.pvd`.

---

## Mesh & Boundary Information

`mesh5_gmsh_p3.m` expects the Gmsh MATLAB structs to contain:

* `TRIANGLES10 (ne x 10)` for velocity; `TRIANGLES6 (ne x 6)` for pressure.
* Boundary lines:

  * Velocity: `LINES4 (nb x 5)` → `[A M1 M2 B flag]`
  * Pressure: `LINES3 (nb x 4)` → `[A M B flag]`

It builds:

* `nodeInfo.velocity.{x,y}`, `nodeInfo.pressure.{x,y}`
* `elemInfo.velElements (ne x 10)`, `elemInfo.presElements (ne x 6)`
* `boundaryInfo.velLine4Elements.flag_<id> = [A M1 M2 B]`
* `boundaryInfo.presLine3Elements.flag_<id> = [A M B]`
* `boundaryInfo.allVelNodes` = union of all velocity boundary nodes

> If your `LINES4` are `[A B M1 M2 flag]`, adapt the reader once; everything downstream will follow.

---

## Solver Details (what happens inside)

* **Time scheme:** First `nBDF1Steps` use BDF1, then BDF2. AB2 predictor defines $\hat u$ for the preconditioner’s convection.
* **Residual:** $R(U) = M U + \alpha\Delta t\,[K_{\text{lin}}U + F_{\text{conv}}(U)] - \text{RHS}$.
* **Matrix-free J·v:** $Mv + \alpha\Delta t\,[K_{\text{lin}}v + J_{\text{conv}}(U)\,v]$ (no convective matrix assembled in the multiply).
* **Preconditioner:** Assemble once per time step with $\hat u$; optionally refreshed mid-Newton if GMRES stalls.
* **BCs:** Enforced strongly on residual rows (inlet profile on $u_x$, no-slip on walls, zero $u_y$ where appropriate). One pressure DOF is pinned.

**Tuning knobs** (in `main1_time_nr_splitfast_gmres.m`):

* `gmres_restart=50`, `gmres_maxit=200`
* Inexact Newton clamp: `eta_min=1e-10`, `eta_max=1e-2`
* Convective Jacobian refresh: `lag_every`, `stall_factor`
* Preconditioner refresh thresholds: `gmres_iter_thresh_total`, `gmres_relres_stall`

---

## Output

* **Snapshots:** `solution_####.vtp` (velocity magnitude), indexed by `solutions.pvd`.
* **Forces (optional):** `force_series.mat` with time / Cd / Cl (when `force_int` is enabled in the driver).
* **On-screen logs:** Newton progress, GMRES counts, crude Cd/Cl preview.

---

## Performance Tips

* **GMRES faster than backslash:** relies on including **advection in the preconditioner**. This repo does that by default: $K_{\text{lin}} + K_{\text{conv}}(\hat u)$.
* **Avoid hot `ismember` paths:** The provided `force_int.m` is written to work with a precomputed **edge lookup** (`boundaryInfo.edgeMap.flag_<id>`) so it does *not* scan elements. If you need forces:

  1. Build the lookup once after meshing (helper provided in comments of `force_int.m`).
  2. Use the new `force_int` (no global searches).
* **Time step:** If Newton iterations grow, try smaller `dt` (or keep more steps in BDF1), or relax the inexact Newton upper tol slightly (`eta_max`).
* **Reordering:** MATLAB’s `lu(A,'vector')` already uses a fill-reducing ordering; leave it on.

---

## Enabling Drag/Lift (optional)

* Keep `force_multiplier` modest (e.g., every 10–50 steps) to minimize cost.
* Use the provided `force_int.m` (edge-mapped version). After calling the mesh reader, add:

  ```matlab
  % after: [nodeInfo, elemInfo, boundaryInfo] = mesh5_gmsh_p3(...)
  boundaryInfo = build_p3_edge_lookup(elemInfo, boundaryInfo); % helper noted in force_int.m
  ```
* In `main1_time_nr_splitfast_gmres.m`, the driver already collects `Cd_series` / `Cl_series` and saves `force_series.mat`.

---

## Known Gotchas

* **Wrong line ordering:** If your `LINES4` endpoint columns are `[A B M1 M2]`, update the reader (`mesh5_gmsh_p3.m`) or the two lines inside `force_int.m` that read endpoints.
* **Pressure pin:** `cfg.corner` must be a valid **P2** index (1…#pressure nodes).
* **NaNs/divergence:** Reduce `dt`, increase BDF1 warm-up (`timeScheme='BDF1(xN)'`), or add small grad-div (`gamma`).

---

## License

MIT — do whatever, credit appreciated.

---

## Citing / Acknowledgement

If this helped your work, a citation or a star on GitHub is always nice ✨.
