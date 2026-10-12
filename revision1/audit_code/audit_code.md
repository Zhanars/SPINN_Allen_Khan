# Audit of the numerical reference and baselines (Symmetry symmetry-4596052, round 1 revision)

Working directory: `D:\Videos\KazNPU\Symmetry_PINN_2026\revision1\audit_code`. Nothing outside this directory was modified (repo code was only imported/read).
All numbers below are produced by the scripts listed in Section 1 (JSON outputs next to them); the tables in Sections 3-6 are generated from those JSON files by `make_tables.py`, not typed by hand.

Conventions. `L2_200` = the repo's discrete norm `sqrt(sum(e^2)*dx)` on the uniform 200-point grid `linspace(-1,1,200)`, `dx=2/199` (what the PINN evaluation uses). `L2_cont` = `sqrt(int_{-1}^{1} e^2 dx)` on a 20001-point grid (trapezoid), with the same interpolation rule as the corresponding method (FD data: linear `np.interp` exactly as in the repo; Chebyshev data: barycentric polynomial interpolation). `Linf_200` = max |e| on the 200-point grid. "Nodal" = compared at the method's own nodes, no interpolation. "Trusted" = the spectral reference of Section 3 (`refsol.reference`).

## 0. Summary

1. **(a) Confirmed.** For `u0(x)=x^2 cos(pi x)`: `u0'(x)=2x cos(pi x) - pi x^2 sin(pi x)`, so `u0'(-1)=+2`, `u0'(+1)=-2` exactly (sympy; finite differences give 2.0039 at h=1e-3 and 2.00004 at h=1e-5). It violates `u_x=0`. The compatible IC `u0c=(cos(pi x)-cos^2(pi x))/2` has `u0c'(+-1)=0` exactly, `u0c(+-1)=-1`, `u0c(0)=u0c(+-0.5)=0`, `u0c''(+-1)=3pi^2/2`.
2. **Trusted reference built and certified** (Chebyshev collocation, Neumann by boundary elimination = index-1 DAE, Radau rtol 1e-11, half-interval even-symmetry reduction). The default settings agree with N-doubling, tolerance tightening, Radau/BDF/LSODA and half- vs full-domain variants to <= 1.2e-10 (T<=1, eps in {0.01,0.05,0.1}, both ICs); an **independent** second-order ghost-node FD solver (Radau, Nx=1000..8000, Richardson-extrapolated) matches it to <= 5.4e-9 nodally (worst case orig IC, eps=0.01, t=0.1; the other 17 cases <= 1e-10). For eps=0.01, T=10 (interfaces formed) N=1024 matches N=1536 to 1.3e-11 and BDF to 2e-9 (N=256 is wrong by 1e-2 there). Certified accuracy: **<= 1e-8 (spectral convergence itself ~1e-11)**. Cost: ~0.2 s for T=1, ~10 s for T=10, eps=0.01.
3. **(b) The current FD reference is substantially wrong near the walls.** Against the trusted solution (orig IC, Nx=500, dt=1e-4, 200-pt grid, linear interpolation as in the repo): Linf error 2.7e-3 / 3.4e-3 / 3.7e-3 at t=0.1 for eps=0.01/0.05/0.1, 1.4e-3 / 1.9e-3 / 2.5e-3 at t=0.5, 5.6e-4 / 8.9e-4 / 1.4e-3 at t=1; L2_cont 1.4e-4 ... 7.3e-4 (t=0.1), 6.3e-5 ... 4.8e-4 (t=1). The error is first order in dx (it drops by 3.4-3.9x when Nx is multiplied by 4, and does not change when dt is divided by 4) and is caused by the boundary rows (`common.py:120-121`). The ghost-node correction reduces it to Linf <= 1.3e-4 (except 7.2e-4 at eps=0.01, t=0.1, orig IC, the under-resolved wall layer), i.e. to about the linear-interpolation floor (L2_200 1.3e-5 .. 3.3e-5, Linf <= 1.2e-4).
4. **The Chebyshev-baseline row (Table 1) is wrong in both directions** because it is scored against the faulty FD reference with a mismatched norm (`l2_error(u_cheb_nodes, interp(FD), dx_ref)` sums over 81 non-uniform nodes with the FD spacing). Stored (N,dt)=(40,2e-3),(60,1e-3),(80,5e-4): L2 = 7.6e-5, 7.8e-5, 8.5e-5; Linf = 5.50e-4, 5.53e-4, 5.56e-4 ("no improvement with N and dt" - an artifact: the Linf is the FD reference's own error, 5.6e-4 at t=1). **Corrected** (vs trusted, proper norms): L2 = 1.69e-4, 8.43e-5, 4.22e-5; Linf = 2.63e-4, 1.30e-4, 6.50e-5, i.e. the baseline error is clean first order in dt (4.2e-6 for N=80, dt=5e-5, 0.46 s). The stored L2 is 2.2x too low for N=40 and 2.0x too high for N=80; the stored Linf is 8.6x too high for N=80. The eps-sweep Chebyshev rows (stored Linf 4.3e-4 .. 3.9e-5) are exactly the FD reference's own error (4.27e-4 .. 3.93e-5); corrected Chebyshev errors there: L2 8.4e-6, 4.2e-6, 1.7e-6, 9.5e-7.
5. **Effect on the PINN numbers is small.** Recomputing the stored seed-42 predictions against the trusted reference changes the Spectral-PINN L2 by +2.2%, -1.1%, -0.1%, +0.5% (t=0.2,0.5,0.8,1.0) and Linf by <= 0.5%; vanilla PINN by <= 0.3% (L2), -5.2% (Linf at t=0.2). Conclusions "Spectral-PINN ~1.4e-3, vanilla ~6.8e-2" stand; but the reference must be replaced for the Chebyshev row, the eps-sweep baseline rows and anything at the 1e-4 level.
6. **(c) The submitted benchmark is nearly eps-independent and nearly diffusion-free on [0,1].** For the orig IC, eps=0.01: max u = 0.056 -> 0.149 (t: 0 -> 1), max|u_x| = 2.5 -> 3.6 (a formed interface would have 70.7), no interface forms until t ~ 3-5 (max|u_x| 15.5 at t=3, 56 at t=5, 70.7 at t=10); max|u(.,1;0.01)-u(.,1;0.02)| = 5.1e-3, max|u(.,1;0.05)-u(.,1;0.1)| = 9.9e-2 (compat IC: 8.4e-3 and 1.5e-1). **Dropping the diffusion term entirely** (pointwise ODE `u'=u-u^3`, closed form) gives L2_200 = 1.0e-3, 9.9e-4, 1.05e-3, 1.30e-3 and Linf = 6.9e-3, 6.1e-3, 4.3e-3, 3.3e-3 at t=0.2,0.5,0.8,1 (eps=0.01, orig IC) - the same size as the Spectral-PINN (L2 2.2e-3, 1.0e-3, 2.5e-3, 1.5e-3; Linf 4.4e-3, 2.2e-3, 5.8e-3, 3.8e-3) and 25-50x better than the vanilla PINN (L2 2.7e-2 ... 6.8e-2). u(.,1;eps=0.01) differs from u(.,1;eps=0.001) by only 3.0e-3 (Linf): the eps-sweep tests essentially one function four times. The trusted solution satisfies the maximum bound principle |u|<=1 in all cases (t>0).
7. **Code issues most likely to bias/invalidate reported results** (Section 7, with file:line): (A3) the BC loss is scaled by 1/K^2 per residual, i.e. 1/K^4 in the loss (effective weight 1e-6 at K=32, 2e-10 at K=256) while the IC fit has wall slope -2: BC effectively unenforced (wall slopes of the trained K=32 network 0.13-0.38 vs 0); (A4) the IC least-squares fit and the PDE collocation use 300 **equispaced** x-points for every K: condition number 41 (K=64), 9e4 (K=96), 1.4e10 (K=128), 4e15 (K>=192); for K=192/256 the fitted IC is wrong by 1.1-1.2 (sup norm) between grid points (wall slope -5e3/-1.8e4), which plausibly explains the K=256 failure (L2 1.5) in the "optimization bottleneck" stress test (not tested by retraining); (A5) the energy penalty is inactive: logged energy loss is nonzero in 1 of 120 sampled epochs (7.6e-8 at epoch 0), with/without penalty differ by <1e-3 relative, and the "zero violations" test cannot fail at T=1 because the true dissipation per sample (1.27e-3-1.34e-3) exceeds twice the largest energy error of the network (4.1e-4); (B3) a no-training oracle test shows the PDE-residual estimator of eps is not intrinsically +40% biased; the sweep's pooled ratio eps_hat/eps_true = 1.394 +- 0.028 (K>=32, n=15) ~ sqrt(2), (eps_hat/eps)^2 = 1.94 ~ 2, constant over a 20x range of eps - a systematic factor, not noise or truncation.
8. **Deliverables**: `refsol.py` (self-test passes in 6-22 s depending on machine load; limit 60 s), this report, scripts + JSON.

## 1. Files in this directory

| file | purpose | output |
|---|---|---|
| `refsol.py` | reusable trusted reference module + self-test (`python refsol.py`) | - |
| `cert_reference.py`, `long_cert.py`, `long_cert2.py`, `independent_check.py` | certification (Section 3) | `cert_reference.json`, `long_cert.json`, `long_cert2.json`, `independent_check.json` |
| `verify_replicas.py` | shows that my copies of the repo FD / Chebyshev solvers are bit-identical (max diff 0.0) to the repo functions and reproduce `data/chebyshev_baseline.json` (L2 8.51097e-5, Linf 5.5627e-4) | stdout |
| `compare_baselines.py` | Section 4 (errors of (i),(ii),(iv) vs trusted) | `compare_baselines.json` |
| `extra_analysis.py {sens,long,energy,epssweep}` | identifiability, long-time FD errors, K-floors, energy, eps-sweep reference errors, corrected Chebyshev rows | `extra_analysis_*.json` |
| `pinn_vs_trusted.py` | stored PINN predictions re-scored vs trusted; K-representation floors | `pinn_vs_trusted.json` |
| `behavior.py`, `ode_baseline.py` | Section 5 | `behavior.json`, `ode_baseline.json`, `ode_vs_pinn.json` |
| `code_checks.py`, `code_checks2.py` | Section 7 evidence | `code_checks.json`, `code_checks2.json` |
| `make_tables.py` | builds the markdown tables from the JSON files | `tables_fragments.json` |

All runs used <= 2 BLAS threads (`OPENBLAS_NUM_THREADS=2`).

## 2. (a) The initial condition does not satisfy the Neumann condition

* `u0(x)=x^2 cos(pi x)`: sympy gives `u0'(x) = 2x cos(pi x) - pi x^2 sin(pi x)`, `u0'(-1)=2`, `u0'(1)=-2`, `u0'(0)=0`. Numerically `(u0(-1+h)-u0(-1))/h = 2.003925 (h=1e-3), 2.0000393 (h=1e-5)`. `u0(+-1)=-1`, `max u0 = 0.05570` at `x=+-0.3428`.
* Consequence for t>0 (BC imposed): a wall layer of thickness ~ eps*sqrt(t) is created. For eps=0.01 the wall value relaxes from -1 to -0.9941 by t=0.1 (eps=0.1: -0.9346); the layer is thinner than the K=32 Chebyshev spacing at the wall (4.8e-3) for t < 0.23 (eps=0.01) and than the FD spacing 4e-3 for t < 0.16.
* `u0c(x)=(cos(pi x)-cos(pi x)^2)/2` (sympy: identical to `sin^2(pi x/2) cos(pi x)`): `u0c'(+-1)=0`, `u0c(+-1)=-1`, `u0c(0)=u0c(+-0.5)=0`, `max u0c = 0.125` at `x=+-1/3`, `u0c''(+-1)=3 pi^2/2`.

## 3. Trusted reference: construction and certification

**Method** (`refsol.reference`). Chebyshev collocation on `x_j=cos(j pi/N)` (same node order as `code/chebyshev_spectral.py`; Weideman-Reddy sine formula for node differences). The Neumann rows `(Du)_0=(Du)_N=0` are used to eliminate the two boundary unknowns (identical to the index-1 DAE with algebraic boundary rows; for t>0 the BC is exact, also for the incompatible IC, whose mismatch at the wall node is confined to a layer of width ~5/N^2). Both ICs are even and the PDE+BC commute with `x -> -x`, so the default solver works on the even half (N/2 unknowns; verified against the full-interval solver, max diff <= 1.8e-13, in the table below). Time integration: `scipy.integrate.solve_ivp`, Radau, analytic Jacobian, rtol 1e-11, atol 1e-13, restarted at each output time (step endpoints, no dense-output error). Evaluation at arbitrary x: barycentric interpolation.

**Default resolution** (`refsol._default_N`): T<=1: N=128 (eps>=0.05), 192 (>=0.02), 256 (>=0.01), 512 (smaller); T>1: N = min(2048, max(256, 2*ceil(5.12/eps))) (= 1024 for eps=0.01).

### 3.1 T=1, eps in {0.1,0.05,0.01}, both ICs (t = 0.1, 0.5, 1; max over t of the Linf difference on a 2001-point grid)

Columns N=64..384: difference to the N=768 solution (same rtol; the N-independent floor of ~3.5e-11 at eps=0.1 (smaller at smaller eps) is, most likely, rounding in the N=768 differentiation matrices rather than an error of the smaller N: the same floor appears between N=96..384 and N=768 alike); `rtol` columns: N=384 with the given rtol vs rtol 1e-12; BDF/LSODA: N=384, rtol 1e-11 vs Radau rtol 1e-12; last column: half-interval vs full-interval solver at N=128.

| IC | eps | N=64 | N=96 | N=128 | N=192 | N=256 | N=384 | rtol 1e-8 | rtol 1e-10 | rtol 1e-11 (default) | BDF | LSODA | half vs full |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| orig | 0.1 | 1.0e-10 | 3.6e-11 | 3.6e-11 | 3.7e-11 | 3.8e-11 | 3.3e-11 | 4.4e-11 | 1.4e-13 | 3.2e-13 | 6.6e-11 | 3.3e-11 | 1.7e-13 |
| orig | 0.05 | 4.2e-08 | 1.8e-11 | 1.8e-11 | 1.8e-11 | 1.9e-11 | 1.6e-11 | 6.1e-11 | 2.9e-13 | 4.0e-14 | 6.7e-11 | 1.7e-11 | 9.1e-14 |
| orig | 0.01 | 1.7e-05 | 6.4e-07 | 5.2e-08 | 3.8e-11 | 3.8e-12 | 3.3e-12 | 1.4e-10 | 6.5e-13 | 3.4e-14 | 1.2e-10 | 5.7e-11 | 1.1e-14 |
| compat | 0.1 | 6.8e-11 | 3.5e-11 | 3.5e-11 | 3.5e-11 | 3.7e-11 | 3.1e-11 | 9.1e-12 | 1.6e-13 | 1.4e-13 | 7.6e-11 | 1.5e-11 | 1.7e-13 |
| compat | 0.05 | 3.0e-08 | 2.7e-11 | 1.8e-11 | 1.8e-11 | 1.9e-11 | 1.6e-11 | 6.6e-12 | 1.0e-13 | 4.2e-14 | 6.3e-11 | 2.6e-11 | 8.6e-14 |
| compat | 0.01 | 9.0e-07 | 1.1e-09 | 3.4e-12 | 3.4e-12 | 3.6e-12 | 3.0e-12 | 2.3e-11 | 7.8e-14 | 1.5e-14 | 9.0e-11 | 5.4e-11 | 1.1e-14 |

Default settings (N=256 at eps=0.01, N=128 at eps>=0.05) are therefore converged to <= 4e-11 in all six cases.

### 3.2 Independent check by a different discretisation

Second-order ghost-node FD (diag `-2`, off-diagonal `2` in the wall rows) in the method of lines, Radau rtol 1e-10 / atol 1e-12 with sparse Jacobian, Nx = 1000, 2000, 4000, 8000, compared **at the common FD nodes** (no interpolation) with the trusted solution (N=256). The error falls by exactly 4x per halving (h^2), and the Richardson extrapolation agrees with the spectral solution to <= 5.4e-9 (worst case orig IC, eps=0.01, t=0.1, where the incompatible IC makes the h^2 asymptotics start later), mostly to 1e-12 - 1e-10. This rules out a common-mode modelling error (BC rows, sign, scaling) in the spectral solver.

| IC | eps | t | FD Nx=1000 | 2000 | 4000 | 8000 | Richardson(2000,4000) | Richardson(4000,8000) |
|---|---|---|---|---|---|---|---|---|
| orig | 0.1 | 0.1 | 1.5e-05 | 3.8e-06 | 9.5e-07 | 2.4e-07 | 8.4e-11 | 4.2e-12 |
| orig | 0.1 | 0.5 | 3.7e-06 | 9.3e-07 | 2.3e-07 | 5.8e-08 | 3.2e-12 | 1.4e-12 |
| orig | 0.1 | 1.0 | 1.8e-06 | 4.4e-07 | 1.1e-07 | 2.8e-08 | 1.6e-12 | 1.7e-12 |
| orig | 0.05 | 0.1 | 3.0e-05 | 7.4e-06 | 1.9e-06 | 4.6e-07 | 6.8e-10 | 4.1e-11 |
| orig | 0.05 | 0.5 | 6.4e-06 | 1.6e-06 | 4.0e-07 | 1.0e-07 | 2.5e-11 | 1.2e-12 |
| orig | 0.05 | 1.0 | 1.7e-06 | 4.2e-07 | 1.0e-07 | 2.6e-08 | 3.2e-12 | 7.3e-13 |
| orig | 0.01 | 0.1 | 1.5e-04 | 3.7e-05 | 9.2e-06 | 2.3e-06 | 8.9e-08 | 5.4e-09 |
| orig | 0.01 | 0.5 | 3.0e-05 | 7.5e-06 | 1.9e-06 | 4.7e-07 | 3.4e-09 | 2.1e-10 |
| orig | 0.01 | 1.0 | 7.9e-06 | 2.0e-06 | 4.9e-07 | 1.2e-07 | 4.3e-10 | 2.7e-11 |
| compat | 0.1 | 0.1 | 1.5e-07 | 3.9e-08 | 9.6e-09 | 2.4e-09 | 8.8e-13 | 8.8e-13 |
| compat | 0.1 | 0.5 | 1.2e-06 | 3.1e-07 | 7.7e-08 | 1.9e-08 | 1.5e-12 | 1.5e-12 |
| compat | 0.1 | 1.0 | 3.3e-06 | 8.2e-07 | 2.1e-07 | 5.1e-08 | 1.6e-12 | 1.6e-12 |
| compat | 0.05 | 0.1 | 4.0e-08 | 1.0e-08 | 2.5e-09 | 6.3e-10 | 4.3e-13 | 4.4e-13 |
| compat | 0.05 | 0.5 | 4.4e-07 | 1.1e-07 | 2.8e-08 | 6.9e-09 | 7.8e-13 | 7.9e-13 |
| compat | 0.05 | 1.0 | 1.9e-06 | 4.8e-07 | 1.2e-07 | 3.0e-08 | 8.3e-13 | 8.3e-13 |
| compat | 0.01 | 0.1 | 1.6e-09 | 4.1e-10 | 1.0e-10 | 2.5e-11 | 9.0e-14 | 9.1e-14 |
| compat | 0.01 | 0.5 | 2.1e-08 | 5.2e-09 | 1.3e-09 | 3.3e-10 | 1.6e-13 | 1.6e-13 |
| compat | 0.01 | 1.0 | 1.3e-07 | 3.1e-08 | 7.9e-09 | 2.0e-09 | 1.8e-13 | 1.8e-13 |

### 3.3 eps=0.01, long times (T=10)

Interfaces of width ~sqrt(2) eps form (the nearest singularity of tanh(x/(sqrt(2) eps)) is at distance 0.022 from the real axis, so N ~ 10/eps is needed for 1e-10). N=256 (the T=1 default) has an error of 1e-3 at t=5 and 1e-2 at t=10; N=1024 is converged.

| IC | comparison (eps=0.01) | t=1 | t=2 | t=3 | t=5 | t=7 | t=10 |
|---|---|---|---|---|---|---|---|
| orig | N=256 vs N=1024 | 4.3e-12 | 4.6e-12 | 1.3e-09 | 1.0e-03 | 4.8e-03 | 9.8e-03 |
| orig | N=384 vs N=1024 | 3.7e-12 | 4.0e-12 | 3.6e-12 | 1.6e-05 | 8.3e-05 | 9.3e-05 |
| orig | N=512 vs N=1024 | 4.1e-12 | 4.4e-12 | 4.0e-12 | 5.0e-07 | 4.3e-06 | 4.9e-06 |
| orig | N=768 vs N=1024 | 6.7e-13 | 8.9e-13 | 3.3e-13 | 2.8e-10 | 6.8e-09 | 8.2e-09 |
| orig | N=1024 vs N=1536 (rtol 1e-12) | 1.1e-11 | 1.2e-11 | 1.2e-11 | 1.1e-11 | 8.1e-12 | 7.5e-12 |
| orig | N=1280 vs N=1536 | 1.2e-11 | 1.3e-11 | 1.3e-11 | 1.2e-11 | 1.2e-11 | 1.1e-11 |
| orig | N=1024, rtol 1e-9 vs N=1536 | 1.2e-11 | 1.1e-11 | 1.2e-11 | 1.1e-11 | 8.0e-12 | 9.6e-12 |
| orig | N=1024, BDF vs N=1536 | 7.8e-10 | 5.9e-10 | 5.4e-10 | 1.4e-09 | 7.3e-10 | 9.8e-10 |
| compat | N=256 vs N=1024 | 4.2e-12 | 5.8e-12 | 3.2e-06 | 2.1e-03 | 3.7e-03 | 5.0e-03 |
| compat | N=384 vs N=1024 | 3.7e-12 | 4.2e-12 | 7.6e-09 | 7.1e-05 | 1.3e-04 | 1.8e-04 |
| compat | N=512 vs N=1024 | 4.1e-12 | 4.6e-12 | 1.0e-11 | 2.6e-06 | 6.0e-06 | 6.3e-06 |
| compat | N=768 vs N=1024 | 3.9e-13 | 9.4e-13 | 6.9e-13 | 1.7e-09 | 5.4e-09 | 5.8e-09 |
| compat | N=1024 vs N=1536 (rtol 1e-12) | 1.1e-11 | 1.2e-11 | 1.2e-11 | 1.1e-11 | 1.0e-11 | 1.0e-11 |
| compat | N=1280 vs N=1536 | 1.2e-11 | 1.3e-11 | 1.3e-11 | 1.2e-11 | 1.3e-11 | 1.3e-11 |
| compat | N=1024, rtol 1e-9 vs N=1536 | 1.1e-11 | 1.2e-11 | 1.2e-11 | 1.1e-11 | 9.7e-12 | 1.3e-11 |
| compat | N=1024, BDF vs N=1536 | 7.6e-10 | 4.7e-10 | 6.2e-10 | 1.9e-09 | 1.6e-09 | 1.2e-09 |

(The 1e-11 level in the N=1536 rows is again the rounding floor of the large matrices; BDF at rtol 1e-10 is less accurate than Radau by construction.)

### 3.4 Other eps used in the repo's eps-sweep (T=1, orig IC)

N=512 vs N=1536: 3.5e-12 (eps=0.005), 1.4e-12 (0.002), 7.0e-13 (0.001); N=256 vs N=384 at eps=0.01: 5.4e-13.

**Certified accuracy** (statement of what was demonstrated; not a rigorous a-posteriori bound): Linf error of `refsol.reference` <= 1e-8 for all of: T<=1, eps in {0.01,0.05,0.1}, both ICs (spectral convergence ~1e-11, independent FD check <= 5.4e-9), eps in {0.005,0.002,0.001} (T=1, orig IC), and ~1e-10 for eps=0.01, T<=10. The module also provides `reference_certified` (N vs 1.5N difference as an error estimate) for other parameter sets.

## 4. (b) Errors of the baselines relative to the trusted reference

Replicas of the repo solvers (with the IC as a parameter) were verified bit-identical to `common.reference_solution_fd` and `chebyshev_spectral.solve_chebyshev` for the orig IC (`verify_replicas.py`: max diff 0.0). (i) = current FD reference (Nx=500, dt=1e-4, rows `1+lam`, `-lam`); (ii) = ghost-node FD (same Nx, dt; wall rows `1+2lam`, `-2lam`); (iv) = `solve_chebyshev(N=80, dt=5e-4)`.

### 4.1 Main table (uniform 200-point grid; FD data linearly interpolated as in the repo, Chebyshev data barycentric)

| IC | eps | t | (i) L2_200 | (i) L2_cont | (i) Linf | (ii) L2_200 | (ii) L2_cont | (ii) Linf | (iv) L2_200 | (iv) L2_cont | (iv) Linf |
|---|---|---|---|---|---|---|---|---|---|---|---|
| orig | 0.01 | 0.1 | 3.80e-04 | 1.35e-04 | 2.68e-03 | 1.03e-04 | 2.62e-05 | 7.18e-04 | 4.14e-06 | 4.12e-06 | 7.62e-06 |
| orig | 0.01 | 0.5 | 2.07e-04 | 1.22e-04 | 1.42e-03 | 2.23e-05 | 1.46e-05 | 1.24e-04 | 2.02e-05 | 2.02e-05 | 2.92e-05 |
| orig | 0.01 | 1.0 | 8.89e-05 | 6.25e-05 | 5.60e-04 | 2.30e-05 | 2.23e-05 | 3.34e-05 | 4.22e-05 | 4.22e-05 | 6.50e-05 |
| orig | 0.05 | 0.1 | 5.92e-04 | 4.63e-04 | 3.41e-03 | 2.22e-05 | 1.42e-05 | 1.24e-04 | 5.25e-06 | 4.82e-06 | 2.06e-05 |
| orig | 0.05 | 0.5 | 4.42e-04 | 3.94e-04 | 1.92e-03 | 1.38e-05 | 1.27e-05 | 2.88e-05 | 1.98e-05 | 1.97e-05 | 3.05e-05 |
| orig | 0.05 | 1.0 | 2.33e-04 | 2.14e-04 | 8.87e-04 | 1.97e-05 | 1.93e-05 | 2.78e-05 | 4.03e-05 | 4.03e-05 | 6.35e-05 |
| orig | 0.1 | 0.1 | 8.27e-04 | 7.29e-04 | 3.72e-03 | 1.65e-05 | 1.25e-05 | 6.90e-05 | 1.17e-05 | 1.10e-05 | 4.08e-05 |
| orig | 0.1 | 0.5 | 7.62e-04 | 7.19e-04 | 2.47e-03 | 1.21e-05 | 1.12e-05 | 1.94e-05 | 2.56e-05 | 2.55e-05 | 3.83e-05 |
| orig | 0.1 | 1.0 | 5.06e-04 | 4.84e-04 | 1.44e-03 | 1.55e-05 | 1.55e-05 | 2.27e-05 | 4.70e-05 | 4.70e-05 | 7.14e-05 |
| compat | 0.01 | 0.1 | 1.76e-05 | 1.62e-05 | 4.73e-05 | 1.63e-05 | 1.62e-05 | 2.42e-05 | 4.46e-06 | 4.46e-06 | 6.47e-06 |
| compat | 0.01 | 0.5 | 2.40e-05 | 2.23e-05 | 6.70e-05 | 2.21e-05 | 2.19e-05 | 3.23e-05 | 2.33e-05 | 2.33e-05 | 2.92e-05 |
| compat | 0.01 | 1.0 | 3.92e-05 | 3.89e-05 | 6.23e-05 | 3.88e-05 | 3.87e-05 | 6.23e-05 | 5.37e-05 | 5.37e-05 | 6.48e-05 |
| compat | 0.05 | 0.1 | 6.34e-05 | 4.65e-05 | 3.85e-04 | 1.61e-05 | 1.60e-05 | 2.40e-05 | 4.35e-06 | 4.35e-06 | 5.95e-06 |
| compat | 0.05 | 0.5 | 9.02e-05 | 7.81e-05 | 4.23e-04 | 2.01e-05 | 1.99e-05 | 2.85e-05 | 2.25e-05 | 2.25e-05 | 3.10e-05 |
| compat | 0.05 | 1.0 | 6.45e-05 | 5.94e-05 | 2.39e-04 | 3.07e-05 | 3.07e-05 | 4.49e-05 | 4.99e-05 | 4.99e-05 | 6.15e-05 |
| compat | 0.1 | 0.1 | 1.63e-04 | 1.38e-04 | 8.10e-04 | 1.54e-05 | 1.53e-05 | 2.33e-05 | 5.37e-06 | 5.37e-06 | 8.91e-06 |
| compat | 0.1 | 0.5 | 2.55e-04 | 2.36e-04 | 9.20e-04 | 1.52e-05 | 1.52e-05 | 1.93e-05 | 2.53e-05 | 2.53e-05 | 3.92e-05 |
| compat | 0.1 | 1.0 | 1.93e-04 | 1.83e-04 | 6.13e-04 | 1.89e-05 | 1.96e-05 | 2.74e-05 | 5.10e-05 | 5.10e-05 | 6.95e-05 |

Reading guide. (i) is wrong at the 6e-5 .. 7e-4 level in L2_cont and 5.6e-4 .. 3.7e-3 in Linf (orig IC; compat IC at eps=0.01 is nearly insensitive, 6e-5). (ii) is at the level of the linear-interpolation floor. (iv) is 4e-6 .. 5.4e-5 (L2) and grows with t (first-order backward Euler, accumulating). `L2_200` overstates `L2_cont` when the error sits in the wall layer (e.g. (i) at t=0.1, eps=0.01: 3.8e-4 vs 1.35e-4: the endpoints of the 200-point grid carry full weight `dx`); for smooth errors both agree to <1%.

### 4.2 Where does the FD error come from? (nodal Linf, so no interpolation involved)

"lin-interp floor": error of linearly interpolating the *exact* solution from the 501 FD nodes to the 200-point grid (the repo's `np.interp`), i.e. the part of any FD-based error that is not the scheme's.

| IC | eps | t | lin-interp floor L2_200 | lin-interp floor Linf | (i) nodal Linf Nx500,dt1e-4 | (i) dt 2.5e-5 | (i) Nx 2000 | (ii) nodal Linf Nx500,dt1e-4 | (ii) dt 2.5e-5 | (ii) Nx 2000 | (i) vs (ii) Linf_200 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| orig | 0.01 | 0.1 | 1.32e-05 | 6.67e-05 | 2.68e-03 | 2.68e-03 | 7.95e-04 | 7.18e-04 | 7.17e-04 | 3.78e-05 | 1.96e-03 |
| orig | 0.01 | 0.5 | 1.78e-05 | 8.29e-05 | 1.42e-03 | 1.42e-03 | 3.83e-04 | 1.24e-04 | 1.23e-04 | 8.19e-06 | 1.30e-03 |
| orig | 0.01 | 1.0 | 2.06e-05 | 3.21e-05 | 5.60e-04 | 5.59e-04 | 1.48e-04 | 3.24e-05 | 3.19e-05 | 1.30e-05 | 5.27e-04 |
| orig | 0.05 | 0.1 | 2.14e-05 | 1.20e-04 | 3.41e-03 | 3.41e-03 | 8.83e-04 | 1.24e-04 | 1.20e-04 | 1.15e-05 | 3.29e-03 |
| orig | 0.05 | 0.5 | 1.52e-05 | 3.49e-05 | 1.92e-03 | 1.91e-03 | 4.91e-04 | 2.88e-05 | 2.64e-05 | 6.13e-06 | 1.89e-03 |
| orig | 0.05 | 1.0 | 1.82e-05 | 2.77e-05 | 8.87e-04 | 8.84e-04 | 2.28e-04 | 1.48e-05 | 7.49e-06 | 1.28e-05 | 8.77e-04 |
| orig | 0.1 | 0.1 | 1.75e-05 | 7.04e-05 | 3.72e-03 | 3.71e-03 | 9.51e-04 | 6.90e-05 | 6.29e-05 | 1.20e-05 | 3.65e-03 |
| orig | 0.1 | 0.5 | 1.34e-05 | 2.38e-05 | 2.47e-03 | 2.46e-03 | 6.27e-04 | 1.94e-05 | 1.60e-05 | 7.77e-06 | 2.45e-03 |
| orig | 0.1 | 1.0 | 1.36e-05 | 1.82e-05 | 1.44e-03 | 1.43e-03 | 3.67e-04 | 1.95e-05 | 9.80e-06 | 1.46e-05 | 1.43e-03 |
| compat | 0.01 | 0.1 | 1.65e-05 | 2.42e-05 | 4.73e-05 | 4.73e-05 | 1.87e-05 | 1.30e-06 | 3.28e-07 | 1.29e-06 | 4.73e-05 |
| compat | 0.01 | 0.5 | 2.12e-05 | 3.10e-05 | 6.70e-05 | 6.69e-05 | 2.04e-05 | 5.86e-06 | 1.48e-06 | 5.84e-06 | 6.69e-05 |
| compat | 0.01 | 1.0 | 3.33e-05 | 5.62e-05 | 3.78e-05 | 3.78e-05 | 1.30e-05 | 1.30e-05 | 3.28e-06 | 1.30e-05 | 3.78e-05 |
| compat | 0.05 | 0.1 | 1.64e-05 | 2.41e-05 | 3.85e-04 | 3.85e-04 | 1.05e-04 | 1.27e-06 | 3.85e-07 | 1.20e-06 | 3.85e-04 |
| compat | 0.05 | 0.5 | 2.00e-05 | 2.87e-05 | 4.23e-04 | 4.23e-04 | 1.10e-04 | 6.77e-06 | 2.61e-06 | 6.24e-06 | 4.23e-04 |
| compat | 0.05 | 1.0 | 2.82e-05 | 4.31e-05 | 2.39e-04 | 2.39e-04 | 6.23e-05 | 1.47e-05 | 8.52e-06 | 1.23e-05 | 2.39e-04 |
| compat | 0.1 | 0.1 | 1.60e-05 | 2.39e-05 | 8.10e-04 | 8.10e-04 | 2.12e-04 | 2.04e-06 | 8.37e-07 | 1.79e-06 | 8.10e-04 |
| compat | 0.1 | 0.5 | 1.68e-05 | 2.29e-05 | 9.20e-04 | 9.19e-04 | 2.36e-04 | 1.04e-05 | 5.93e-06 | 7.94e-06 | 9.19e-04 |
| compat | 0.1 | 1.0 | 1.85e-05 | 2.57e-05 | 6.13e-04 | 6.11e-04 | 1.58e-04 | 2.23e-05 | 1.49e-05 | 1.42e-05 | 6.12e-04 |

* (i): the nodal error is unchanged when dt is divided by 4 and drops by 3.4-3.9x when Nx is multiplied by 4: it is a **first-order spatial error from the wall rows**, not a time error. With the orig IC it is 2.7e-3 .. 3.7e-3 at t=0.1.
* (ii): unchanged by dt at t<=0.5 and 6-19x smaller at Nx=2000 for the orig IC at t=0.1 (eps=0.01: 7.2e-4 -> 3.8e-5; eps=0.1: 6.9e-5 -> 1.2e-5), so the remaining 7.2e-4 is the O(h^2) resolution error of the incompatible-IC wall layer, not a bug; at t=1 what is left at Nx=2000 (1.3e-5 ... 1.5e-5) is the backward-Euler time error (compare the Chebyshev baseline with the same scheme: 6.5e-5 at dt=5e-4).
* The linear interpolation (`np.interp`) costs L2_200 1.3e-5 .. 3.3e-5 and Linf up to 1.2e-4 whatever the FD scheme.

### 4.3 Chebyshev baseline: repo-style evaluation vs corrected

"repo-style" = `l2_error(u_cheb, np.interp(x_cheb, x_FD, u_FD), dx_FD)` and `linf_error` (what `run_experiments.py:86-87`, `chebyshev_spectral.py:84-86` do) for (iv) vs the current FD reference. "corrected" = (iv) vs trusted, L2 by Clenshaw-Curtis quadrature on the nodes (equal to the barycentric `L2_cont`), Linf at the nodes.

| IC | eps | t | repo-style L2 (vs FD ref, dx_ref weights) | repo-style Linf (vs FD ref) | corrected L2 (Clenshaw-Curtis, vs trusted) | corrected L2_cont (barycentric) | corrected Linf (nodes, vs trusted) |
|---|---|---|---|---|---|---|---|
| orig | 0.01 | 0.1 | 3.06e-04 | 2.68e-03 | 4.03e-06 | 4.12e-06 | 6.41e-06 |
| orig | 0.01 | 0.5 | 1.97e-04 | 1.42e-03 | 2.02e-05 | 2.02e-05 | 2.89e-05 |
| orig | 0.01 | 1.0 | 8.51e-05 | 5.56e-04 | 4.22e-05 | 4.22e-05 | 6.42e-05 |
| orig | 0.05 | 0.1 | 5.75e-04 | 3.39e-03 | 4.82e-06 | 4.82e-06 | 2.06e-05 |
| orig | 0.05 | 0.5 | 3.88e-04 | 1.90e-03 | 1.97e-05 | 1.97e-05 | 3.03e-05 |
| orig | 0.05 | 1.0 | 1.90e-04 | 8.71e-04 | 4.03e-05 | 4.03e-05 | 6.34e-05 |
| orig | 0.1 | 0.1 | 7.36e-04 | 3.68e-03 | 1.10e-05 | 1.10e-05 | 4.08e-05 |
| orig | 0.1 | 0.5 | 5.87e-04 | 2.44e-03 | 2.55e-05 | 2.55e-05 | 3.82e-05 |
| orig | 0.1 | 1.0 | 3.62e-04 | 1.41e-03 | 4.70e-05 | 4.70e-05 | 7.09e-05 |
| compat | 0.01 | 0.1 | 8.63e-06 | 4.73e-05 | 4.46e-06 | 4.46e-06 | 6.41e-06 |
| compat | 0.01 | 0.5 | 1.33e-05 | 6.68e-05 | 2.33e-05 | 2.33e-05 | 2.88e-05 |
| compat | 0.01 | 1.0 | 1.62e-05 | 5.51e-05 | 5.37e-05 | 5.37e-05 | 6.44e-05 |
| compat | 0.05 | 0.1 | 6.06e-05 | 3.85e-04 | 4.35e-06 | 4.35e-06 | 5.94e-06 |
| compat | 0.05 | 0.5 | 8.05e-05 | 4.20e-04 | 2.25e-05 | 2.25e-05 | 3.09e-05 |
| compat | 0.05 | 1.0 | 4.98e-05 | 2.35e-04 | 4.99e-05 | 4.99e-05 | 6.07e-05 |
| compat | 0.1 | 0.1 | 1.51e-04 | 8.09e-04 | 5.37e-06 | 5.37e-06 | 8.88e-06 |
| compat | 0.1 | 0.5 | 2.05e-04 | 9.12e-04 | 2.53e-05 | 2.53e-05 | 3.87e-05 |
| compat | 0.1 | 1.0 | 1.45e-04 | 6.02e-04 | 5.10e-05 | 5.10e-05 | 6.96e-05 |

The row orig / eps=0.01 / t=1.0 reproduces the stored `data/chebyshev_baseline.json` row N=80 (8.51e-5, 5.56e-4). Table 1 of the manuscript and the eps-sweep rows, corrected (all vs trusted, orig IC, eps=0.01 unless stated; timings are from my LU-based replica, N=80, dt=5e-4: 0.05 s; dt=5e-5: 0.46 s):

| eps | N | dt | stored L2 (repo) | stored Linf (repo) | corrected L2_cont | corrected Linf_200 | corrected L2_200 |
|---|---|---|---|---|---|---|---|
| 0.01 | 40 | 0.002 | 7.58e-05 | 5.50e-04 | 1.69e-04 | 2.63e-04 | 1.69e-04 |
| 0.01 | 60 | 0.001 | 7.78e-05 | 5.53e-04 | 8.43e-05 | 1.30e-04 | 8.43e-05 |
| 0.01 | 80 | 0.0005 | 8.51e-05 | 5.56e-04 | 4.22e-05 | 6.50e-05 | 4.22e-05 |
| 0.01 | 80 | 5e-05 | - | - | 4.22e-06 | 6.51e-06 | 4.22e-06 |
| 0.01 | 400 | 0.0001 | 1.17e-04 | 4.26e-04 | 8.43e-06 | 1.30e-05 | 8.44e-06 |
| 0.005 | 400 | 5e-05 | 3.39e-05 | 2.04e-04 | 4.22e-06 | 6.52e-06 | 4.22e-06 |
| 0.002 | 400 | 2e-05 | 7.05e-06 | 7.93e-05 | 1.72e-06 | 2.82e-06 | 1.72e-06 |
| 0.001 | 400 | 1e-05 | 2.66e-06 | 3.93e-05 | 9.53e-07 | 1.84e-06 | 9.47e-07 |

The stored Linf in the eps-sweep (4.26e-4, 2.04e-4, 7.93e-5, 3.93e-5) equals the FD-reference error vs trusted computed below (4.27e-4, 2.04e-4, 7.93e-5, 3.93e-5), i.e. the "Chebyshev error decreasing with eps" is the FD reference refinement `reference_nx_for_eps` showing through.

| eps | nx_ref (repo) | FD-current ref L2_200 vs trusted | FD-current ref Linf vs trusted |
|---|---|---|---|
| 0.01 | 668 | 6.76e-05 | 4.27e-04 |
| 0.005 | 1334 | 3.05e-05 | 2.04e-04 |
| 0.002 | 3334 | 1.41e-05 | 7.93e-05 |
| 0.001 | 6668 | 1.01e-05 | 3.93e-05 |

### 4.4 Stored seed-42 PINN predictions re-scored against the trusted reference

`data/spectral_pinn.json`, `spectral_pinn_energy.json`, `vanilla_pinn.json` contain `u_pred` on the 200-point grid for t = 0.2, 0.5, 0.8, 1.0 (eps=0.01, orig IC). The stored errors reproduce exactly when recomputed against the repo's FD reference; against the trusted reference:

| run (seed 42) | t | stored L2 (vs FD ref) | stored Linf | L2 vs trusted | Linf vs trusted | dL2 | dLinf |
|---|---|---|---|---|---|---|---|
| spectral_pinn | 0.2 | 2.120e-03 | 4.437e-03 | 2.167e-03 | 4.437e-03 | +2.2% | -0.0% |
| spectral_pinn | 0.5 | 1.037e-03 | 2.217e-03 | 1.026e-03 | 2.221e-03 | -1.1% | +0.2% |
| spectral_pinn | 0.8 | 2.492e-03 | 5.822e-03 | 2.489e-03 | 5.835e-03 | -0.1% | +0.2% |
| spectral_pinn | 1.0 | 1.451e-03 | 3.784e-03 | 1.458e-03 | 3.802e-03 | +0.5% | +0.5% |
| spectral_pinn_energy | 0.2 | 2.120e-03 | 4.437e-03 | 2.168e-03 | 4.437e-03 | +2.2% | -0.0% |
| spectral_pinn_energy | 0.5 | 1.038e-03 | 2.219e-03 | 1.027e-03 | 2.223e-03 | -1.1% | +0.2% |
| spectral_pinn_energy | 0.8 | 2.494e-03 | 5.825e-03 | 2.491e-03 | 5.839e-03 | -0.1% | +0.2% |
| spectral_pinn_energy | 1.0 | 1.452e-03 | 3.787e-03 | 1.459e-03 | 3.804e-03 | +0.5% | +0.5% |
| vanilla_pinn | 0.2 | 2.674e-02 | 4.571e-02 | 2.665e-02 | 4.333e-02 | -0.3% | -5.2% |
| vanilla_pinn | 0.5 | 4.197e-02 | 6.997e-02 | 4.197e-02 | 6.998e-02 | -0.0% | +0.0% |
| vanilla_pinn | 0.8 | 5.658e-02 | 8.397e-02 | 5.660e-02 | 8.399e-02 | +0.0% | +0.0% |
| vanilla_pinn | 1.0 | 6.761e-02 | 8.581e-02 | 6.763e-02 | 8.583e-02 | +0.0% | +0.0% |

The multi-seed runs (`forward_multiseed.json`) store no predictions, so only the seed-42 runs can be re-scored.

Representation floor of the PINN's basis (least-squares degree-K Chebyshev fit on 300 equispaced points of the trusted solution, evaluated on the 200-point grid; L2_200 / Linf):

| K | t=0.1 | t=0.2 | t=0.5 | t=0.8 | t=1.0 |
|---|---|---|---|---|---|
| 16 | 3.9e-04 / 1.9e-03 | 4.0e-04 / 1.9e-03 | 3.4e-04 / 1.1e-03 | 8.4e-04 / 2.1e-03 | 2.0e-03 / 4.6e-03 |
| 24 | 1.6e-04 / 5.6e-04 | 1.5e-04 / 5.8e-04 | 6.7e-05 / 2.4e-04 | 9.5e-05 / 2.0e-04 | 2.4e-04 / 3.3e-04 |
| 32 | 5.6e-05 / 1.8e-04 | 3.9e-05 / 1.4e-04 | 6.5e-06 / 2.5e-05 | 6.6e-06 / 1.7e-05 | 1.6e-05 / 1.9e-05 |
| 48 | 5.3e-06 / 1.9e-05 | 5.4e-06 / 3.4e-05 | 9.7e-07 / 2.4e-06 | 4.7e-07 / 2.6e-06 | 2.4e-07 / 9.4e-07 |
| 64 | 2.2e-06 / 1.5e-05 | 1.0e-06 / 7.0e-06 | 2.1e-07 / 1.4e-06 | 6.8e-08 / 4.8e-07 | 2.6e-08 / 1.7e-07 |

At t=1 the K=32 floor (1.6e-5 / 1.9e-5) is ~90x below the trained network's error (1.46e-3 / 3.8e-3): the accuracy is limited by training, not by K. K=16 would not suffice (2.0e-3 at t=1).

### 4.5 FD reference at long times (eps=0.01, repo parameters Nx=500, dt=1e-4, as in `stress_time_sweep.py`)

| IC | t | (i) L2_200 | (i) Linf | (ii) L2_200 | (ii) Linf |
|---|---|---|---|---|---|
| orig | 1.0 | 8.89e-05 | 5.60e-04 | 2.30e-05 | 3.34e-05 |
| orig | 2.0 | 6.47e-05 | 1.27e-04 | 6.32e-05 | 1.27e-04 |
| orig | 3.0 | 1.60e-04 | 4.76e-04 | 1.60e-04 | 4.76e-04 |
| orig | 5.0 | 5.52e-04 | 2.47e-03 | 5.52e-04 | 2.47e-03 |
| orig | 7.0 | 5.69e-04 | 2.54e-03 | 5.69e-04 | 2.54e-03 |
| orig | 10.0 | 5.66e-04 | 2.48e-03 | 5.66e-04 | 2.48e-03 |
| compat | 1.0 | 3.92e-05 | 6.23e-05 | 3.88e-05 | 6.23e-05 |
| compat | 2.0 | 1.10e-04 | 2.78e-04 | 1.10e-04 | 2.78e-04 |
| compat | 3.0 | 3.02e-04 | 1.17e-03 | 3.02e-04 | 1.17e-03 |
| compat | 5.0 | 5.65e-04 | 2.68e-03 | 5.65e-04 | 2.68e-03 |
| compat | 7.0 | 5.55e-04 | 2.58e-03 | 5.55e-04 | 2.58e-03 |
| compat | 10.0 | 5.53e-04 | 2.55e-03 | 5.53e-04 | 2.55e-03 |

From t >= 2 the wall rows no longer matter (current = ghost) and the error (2.5e-3 Linf, 5.5e-4 L2 at t>=5) is the Nx=500 resolution of the pinned interfaces; negligible against the reported PINN errors of ~1 at t>=5.

K-representation floors at long times (orig IC, eps=0.01; "L2_200 / max|fit|", LS fit on 300 equispaced points as the code does; a value of max|fit|>1 means the series violates |u|<=1):

| K (LS fit on 300 equispaced pts, as the code does) | t=1 | t=2 | t=3 | t=5 | t=7 | t=10 |
|---|---|---|---|---|---|---|
| 32 | 1.56e-05 / 0.997 | 1.53e-03 / 1.000 | 1.93e-02 / 1.016 | 1.35e-01 / 1.168 | 1.56e-01 / 1.183 | 1.57e-01 / 1.180 |
| 64 | 2.56e-08 / 0.997 | 8.00e-06 / 0.999 | 1.51e-03 / 1.003 | 4.49e-02 / 1.044 | 6.23e-02 / 1.097 | 6.28e-02 / 1.098 |
| 128 | 1.94e-08 / 0.997 | 1.39e-02 / 0.999 | 2.77e+02 / 1951.634 | 4.76e+04 / 335862.295 | 6.91e+04 / 487668.407 | 7.49e+04 / 528076.270 |
| 256 | 1.62e-01 / 1.212 | 1.60e-01 / 1.225 | 1.01e+01 / 41.471 | 9.12e+06 / 37023712.451 | 2.81e+07 / 114449995.940 | 2.99e+07 / 121563006.687 |

K=32 cannot do better than L2 0.135-0.157 at t>=5 (Gibbs overshoot 1.17-1.18), K=64 0.045-0.063; the reported PINN errors (0.90, 0.99 at t=5, 10) are ~6x above the K=32 floor. The K>=128 rows are the ill-conditioning of the equispaced fit (Section 7, A4), not properties of the solution.

## 5. (c) Behaviour of the solutions; Neumann-compatible IC; eps-sensitivity

### 5.1 eps = 0.01 up to T = 10 (both ICs)

Interface = sign change of u for x>0 (u is even); `max abs u_x` of a fully formed tanh interface is 1/(sqrt(2) eps) = 70.7; E = Ginzburg-Landau energy (Clenshaw-Curtis; the steady state with two interfaces has E = 2*(2 sqrt2/3) eps = 0.018856).

| IC | t | max u | min u | max abs u_x | E(t) | zero x>0 | u(0,t) |
|---|---|---|---|---|---|---|---|
| orig | 0.0 | 0.0557 | -1.0000 | 2.5 | 0.40044 | 0.5000 | 0.0000 |
| orig | 1.0 | 0.1490 | -0.9967 | 3.6 | 0.33455 | 0.4992 | 0.0005 |
| orig | 2.0 | 0.3772 | -0.9993 | 6.9 | 0.26638 | 0.4984 | 0.0029 |
| orig | 3.0 | 0.7400 | -0.9999 | 15.5 | 0.16485 | 0.4976 | 0.0119 |
| orig | 5.0 | 0.9923 | -1.0000 | 56.1 | 0.05144 | 0.4968 | 0.1444 |
| orig | 7.0 | 0.9999 | -1.0000 | 69.8 | 0.02029 | 0.4968 | 0.8218 |
| orig | 10.0 | 1.0000 | -1.0000 | 70.7 | 0.01886 | 0.4968 | 0.9996 |
| compat | 0.0 | 0.1250 | -1.0000 | 2.8 | 0.36035 | 0.5000 | 0.0000 |
| compat | 1.0 | 0.3223 | -0.9998 | 5.1 | 0.28914 | 0.4994 | 0.0013 |
| compat | 2.0 | 0.6770 | -0.9999 | 11.7 | 0.18771 | 0.4988 | 0.0072 |
| compat | 3.0 | 0.9277 | -1.0000 | 27.7 | 0.09548 | 0.4982 | 0.0294 |
| compat | 5.0 | 0.9985 | -1.0000 | 64.7 | 0.03482 | 0.4979 | 0.3359 |
| compat | 7.0 | 1.0000 | -1.0000 | 70.4 | 0.01895 | 0.4979 | 0.9580 |
| compat | 10.0 | 1.0000 | -1.0000 | 70.7 | 0.01886 | 0.4979 | 0.9999 |

* No interface has formed at T=1 in either IC (max|u_x| 3.6 and 5.1): the region |x|>0.5 with u<0 is driven toward -1 by the reaction (u(0.75,1) = -0.76 orig / -0.90 compat), the small positive hump (0.056 resp. 0.125) grows like e^t (reaction dominates), and the zero crossing stays at x ~ 0.499. Sharp interfaces (max|u_x| > 50) appear around t ~ 5 and the profile is the final one (u=+-1, interface pinned near x = +-0.497 (orig) / +-0.498 (compat), E -> 0.01886 = 2 interfaces) by t ~ 7-10. A meaningful interface test needs T >= 5; the long-time end state is static.
* Max bound principle: `max|u| <= 1` holds for all t>0, all eps, both ICs (`mbp_ok` in `behavior.json`); u=-1 at the walls at t=0 only touches the bound.
* Compat IC vs orig IC: same sign structure, larger hump (0.125 vs 0.056), interface formation ~1 time unit earlier (max|u_x| at t=3: 27.7 vs 15.5), larger eps-sensitivity at t=1.

### 5.2 eps-dependence on [0,1]

| IC | eps | t | max u | min u | max abs u_x | E | zero x>0 |
|---|---|---|---|---|---|---|---|
| orig | 0.01 | 0.0 | 0.0557 | -1.0000 | 2.55 | 0.40044 | 0.5000 |
| orig | 0.01 | 0.5 | 0.0913 | -0.9939 | 2.79 | 0.36705 | 0.4996 |
| orig | 0.01 | 1.0 | 0.1490 | -0.9967 | 3.57 | 0.33455 | 0.4992 |
| orig | 0.02 | 0.0 | 0.0557 | -1.0000 | 2.55 | 0.40109 | 0.5000 |
| orig | 0.02 | 0.5 | 0.0905 | -0.9872 | 2.78 | 0.36773 | 0.4984 |
| orig | 0.02 | 1.0 | 0.1464 | -0.9929 | 3.53 | 0.33548 | 0.4968 |
| orig | 0.05 | 0.0 | 0.0557 | -1.0000 | 2.55 | 0.40567 | 0.5000 |
| orig | 0.05 | 0.5 | 0.0849 | -0.9641 | 2.70 | 0.37226 | 0.4899 |
| orig | 0.05 | 1.0 | 0.1289 | -0.9782 | 3.26 | 0.34115 | 0.4798 |
| orig | 0.1 | 0.0 | 0.0557 | -1.0000 | 2.55 | 0.42201 | 0.5000 |
| orig | 0.1 | 0.5 | 0.0674 | -0.9148 | 2.46 | 0.38593 | 0.4587 |
| orig | 0.1 | 1.0 | 0.0821 | -0.9391 | 2.62 | 0.35356 | 0.4160 |
| compat | 0.01 | 0.0 | 0.1250 | -1.0000 | 2.76 | 0.36035 | 0.5000 |
| compat | 0.01 | 0.5 | 0.2028 | -0.9997 | 3.59 | 0.32628 | 0.4997 |
| compat | 0.01 | 1.0 | 0.3223 | -0.9998 | 5.11 | 0.28914 | 0.4994 |
| compat | 0.02 | 0.0 | 0.1250 | -1.0000 | 2.76 | 0.36109 | 0.5000 |
| compat | 0.02 | 0.5 | 0.2010 | -0.9989 | 3.57 | 0.32746 | 0.4987 |
| compat | 0.02 | 1.0 | 0.3171 | -0.9992 | 5.02 | 0.29137 | 0.4975 |
| compat | 0.05 | 0.0 | 0.1250 | -1.0000 | 2.76 | 0.36628 | 0.5000 |
| compat | 0.05 | 0.5 | 0.1892 | -0.9931 | 3.43 | 0.33516 | 0.4921 |
| compat | 0.05 | 1.0 | 0.2826 | -0.9947 | 4.46 | 0.30462 | 0.4844 |
| compat | 0.1 | 0.0 | 0.1250 | -1.0000 | 2.76 | 0.38478 | 0.5000 |
| compat | 0.1 | 0.5 | 0.1524 | -0.9715 | 3.01 | 0.35681 | 0.4677 |
| compat | 0.1 | 1.0 | 0.1889 | -0.9754 | 3.30 | 0.33064 | 0.4360 |

### 5.3 Sensitivity to eps (Linf / L2 of the difference of trusted solutions)

| IC | eps pair | t=0.1 (Linf / L2) | t=0.25 | t=0.5 | t=1.0 |
|---|---|---|---|---|---|
| orig | 0.01 vs 0.02 | 6.07e-03 / 6.90e-04 | 7.36e-03 / 1.15e-03 | 6.65e-03 / 1.67e-03 | 5.13e-03 / 3.64e-03 |
| orig | 0.05 vs 0.1 | 3.43e-02 / 9.74e-03 | 4.61e-02 / 1.93e-02 | 4.93e-02 / 3.39e-02 | 9.92e-02 / 7.32e-02 |
| orig | 0.01 vs 0.1 | 5.94e-02 / 1.35e-02 | 7.75e-02 / 2.63e-02 | 7.91e-02 / 4.58e-02 | 1.38e-01 / 1.01e-01 |
| compat | 0.01 vs 0.02 | 3.63e-04 / 3.30e-04 | 9.81e-04 / 8.54e-04 | 2.52e-03 / 1.94e-03 | 8.36e-03 / 5.60e-03 |
| compat | 0.05 vs 0.1 | 9.00e-03 / 8.08e-03 | 2.30e-02 / 2.02e-02 | 5.49e-02 / 4.27e-02 | 1.49e-01 / 1.03e-01 |
| compat | 0.01 vs 0.1 | 1.19e-02 / 1.07e-02 | 3.08e-02 / 2.70e-02 | 7.45e-02 / 5.79e-02 | 2.10e-01 / 1.45e-01 |

Asked-for numbers: `max|u(.,1;0.01)-u(.,1;0.02)|` = **5.13e-3 (orig), 8.36e-3 (compat)**; `max|u(.,1;0.05)-u(.,1;0.1)|` = **9.92e-2 (orig), 1.49e-1 (compat)**. Where the difference lives:

| IC | eps pair | t | Linf | argmax x | Linf on |x|<=0.9 | Linf on |x|<=0.5 |
|---|---|---|---|---|---|---|
| orig | 0.01 vs 0.02 | 0.1 | 6.07e-03 | -1.000 | 2.21e-04 | 2.08e-04 |
| orig | 0.01 vs 0.02 | 1.0 | 5.13e-03 | -0.519 | 5.13e-03 | 5.08e-03 |
| orig | 0.05 vs 0.1 | 0.1 | 3.43e-02 | -1.000 | 5.45e-03 | 5.14e-03 |
| orig | 0.05 vs 0.1 | 1.0 | 9.92e-02 | 0.484 | 9.92e-02 | 9.92e-02 |
| compat | 0.01 vs 0.02 | 0.1 | 3.63e-04 | -1.000 | 3.38e-04 | 3.38e-04 |
| compat | 0.01 vs 0.02 | 1.0 | 8.36e-03 | 0.468 | 8.36e-03 | 8.36e-03 |
| compat | 0.05 vs 0.1 | 0.1 | 9.00e-03 | -1.000 | 8.31e-03 | 8.24e-03 |
| compat | 0.05 vs 0.1 | 1.0 | 1.49e-01 | -0.427 | 1.49e-01 | 1.49e-01 |

For the orig IC at t=0.1 the eps-signal sits in the wall layer created by the IC/BC conflict (x = -1; in |x|<=0.9 it is only 2.2e-4); at t=1 it sits near the zero crossing x ~ +-0.5 (5.1e-3). Longer times (eps 0.01 vs 0.02):

| IC (eps 0.01 vs 0.02) | t=1 | t=2 | t=3 | t=5 | t=10 |
|---|---|---|---|---|---|
| orig | 5.13e-03 / 3.64e-03 | 2.67e-02 / 1.51e-02 | 1.02e-01 / 4.23e-02 | 4.44e-01 / 1.45e-01 | 5.41e-01 / 1.23e-01 |
| compat | 8.36e-03 / 5.60e-03 | 4.54e-02 / 2.23e-02 | 1.59e-01 / 5.73e-02 | 4.53e-01 / 1.45e-01 | 4.43e-01 / 1.04e-01 |

### 5.4 The benchmark is (almost) a pointwise ODE: error of dropping the diffusion term

`u(x,t)=u0 e^t / sqrt(1+u0^2 (e^{2t}-1))` solves `u'=u-u^3` (eps=0). Error vs the trusted solution on the 200-point grid (L2_200 / Linf):

| IC | eps | t=0.1 (L2_200 / Linf) | t=0.5 | t=1.0 |
|---|---|---|---|---|
| orig | 0.001 | 8.30e-05 / 5.85e-04 | 8.37e-05 / 5.90e-04 | 4.52e-05 / 3.08e-04 |
| orig | 0.01 | 8.42e-04 / 5.92e-03 | 9.86e-04 / 6.13e-03 | 1.30e-03 / 3.28e-03 |
| orig | 0.05 | 5.48e-03 / 3.11e-02 | 1.32e-02 / 3.59e-02 | 2.92e-02 / 4.08e-02 |
| orig | 0.1 | 1.53e-02 / 6.54e-02 | 4.71e-02 / 8.52e-02 | 1.02e-01 / 1.39e-01 |
| compat | 0.001 | 1.11e-06 / 1.21e-06 | 6.53e-06 / 8.46e-06 | 1.89e-05 / 2.83e-05 |
| compat | 0.01 | 1.11e-04 / 1.21e-04 | 6.52e-04 / 8.45e-04 | 1.89e-03 / 2.82e-03 |
| compat | 0.05 | 2.76e-03 / 3.02e-03 | 1.59e-02 / 2.05e-02 | 4.41e-02 / 6.53e-02 |
| compat | 0.1 | 1.09e-02 / 1.20e-02 | 5.86e-02 / 7.54e-02 | 1.47e-01 / 2.13e-01 |

For eps=0.01 and the orig IC at the PINN evaluation times (L2_200 / Linf; trusted reference): ODE 9.9e-4 / 6.9e-3 (t=0.2), 9.9e-4 / 6.1e-3 (0.5), 1.05e-3 / 4.3e-3 (0.8), 1.30e-3 / 3.3e-3 (1.0); Spectral-PINN seed 42: 2.2e-3 / 4.4e-3, 1.0e-3 / 2.2e-3, 2.5e-3 / 5.8e-3, 1.5e-3 / 3.8e-3; vanilla PINN 2.7e-2, 4.2e-2, 5.7e-2, 6.8e-2 (L2). The ODE error is concentrated at the wall (argmax x=-1). Also `max|u(.,1;eps) - u(.,1;0.001)|` for eps = 0.01, 0.005, 0.002 = 3.0e-3, 1.3e-3, 3.1e-4 (orig) and 2.8e-3, 6.8e-4, 8.5e-5 (compat): the four eps values of `stress_eps_sweep.py` (0.01 ... 0.001) give solutions within 3e-3 of each other at t=1, so constant PINN errors across that sweep (vanilla 6.75e-2, 6.79e-2, 6.81e-2, 6.81e-2; Spectral 1.45e-3, 9.9e-4, 9.8e-4, 9.8e-4) are expected and say nothing about stiffness.

### 5.5 Identifiability of eps from the repo's observation design (40 uniform x-points incl. walls, t = 0.2,...,1.0; 200 values)

RMS change of the observation vector when eps is changed by the given factor (trusted solutions), and RMS error of the observations the repo actually generates (`make_observations`: current FD reference + `np.interp`) relative to the trusted solution:

| IC | eps_true | RMS error of repo's FD-generated data vs trusted | RMS data change for eps x0.95 | x1.05 | x1.25 | x1.4 | x2.0 |
|---|---|---|---|---|---|---|---|
| orig | 0.01 | 3.32e-04 | 8.01e-05 | 8.18e-05 | 4.26e-04 | 7.04e-04 | 1.99e-03 |
| orig | 0.1 | 6.25e-04 | 4.00e-03 | 4.10e-03 | 2.13e-02 | 3.47e-02 | 8.86e-02 |
| compat | 0.01 | (n/a) | 7.59e-05 | 7.97e-05 | 4.37e-04 | 7.45e-04 | 2.32e-03 |
| compat | 0.1 | (n/a) | 5.15e-03 | 5.22e-03 | 2.66e-02 | 4.27e-02 | 1.04e-01 |

At eps=0.01 a 25% (40%) change of eps changes the data by RMS 4.3e-4 (7.0e-4), the same size as the repo's own data-generation error (RMS 3.3e-4, Linf 2.4e-3) and comparable to the forward error of the trained Spectral-PINN (RMS ~1.0e-3 = L2 1.45e-3/sqrt(2)): the inverse problem at eps=0.01 is practically ill-posed with this design. At eps=0.1 the signal is 6x (5% change of eps) to 56x (40% change) above the data error, so the +40% bias reported for eps_true=0.1 cannot be explained by data-generation error or by weak sensitivity.

## 6. Energy diagnostics

Energy of the trained Spectral-PINN (stored `E_traj`, 51 time points, eps=0.01, orig IC) vs the trusted energy (Clenshaw-Curtis):

| run | E_pinn(0) | E_true(0) | E_pinn(1) | E_true(1) | max abs E error over 51 pts | true dissipation per 0.02 step | violations counted |
|---|---|---|---|---|---|---|---|
| spectral_pinn | 0.40046 | 0.40044 | 0.33477 | 0.33455 | 4.13e-04 | 1.274e-03 .. 1.340e-03 | 0 |
| spectral_pinn_energy | 0.40046 | 0.40044 | 0.33477 | 0.33455 | 4.13e-04 | 1.274e-03 .. 1.340e-03 | 0 |

* True energy decreases monotonically by 1.27e-3 .. 1.34e-3 per 0.02-step. A "violation" (E[i+1] > E[i] + 1e-10) would need the network's energy error to change by more than that between consecutive points, but the largest energy error over the whole trajectory is 4.1e-4. The check **cannot fail** for any network of this accuracy at T=1; "zero violations" carries no information there. (At T=5 and T=10 the network's E_last = 0.249 and 0.263 compare with the true E(5)=0.0514 (4.8x) and E(10)=0.0189 (14x); the true energy reaches its plateau at t~7.)
* The two Spectral-PINN runs (with / without the penalty) have E trajectories differing by at most 4.5e-7 and L2 errors by <= 0.1% (Figure 1: identical curves).
* Quadrature: the 300-point uniform trapezoid used for `E` (`spectral_pinn.py:115-120`) is accurate for resolved fields; checked on the trusted solution: relative difference to Clenshaw-Curtis <= 1e-15 for t<=1, 1.3e-8 / 3.8e-7 at t=5 / 10 (orig, eps=0.01), <= 3e-6 (compat). Not a source of bias; the problem at long times is that the K=32 series is not a resolved representation.

| case | E (Clenshaw-Curtis) | E (trapezoid, 300 uniform pts) | rel. diff |
|---|---|---|---|
| orig_eps0.01_t0.5 | 0.36705070 | 0.36705070 | 3.0e-16 |
| orig_eps0.01_t1.0 | 0.33455451 | 0.33455451 | 0.0e+00 |
| orig_eps0.01_t5.0 | 0.05143686 | 0.05143686 | 1.3e-08 |
| orig_eps0.01_t10.0 | 0.01885619 | 0.01885620 | 3.8e-07 |
| orig_eps0.1_t0.5 | 0.38592945 | 0.38592945 | 0.0e+00 |
| orig_eps0.1_t1.0 | 0.35356399 | 0.35356399 | 0.0e+00 |
| compat_eps0.01_t0.5 | 0.32628326 | 0.32628326 | 3.4e-16 |
| compat_eps0.01_t1.0 | 0.28913983 | 0.28913983 | 1.9e-16 |
| compat_eps0.01_t5.0 | 0.03482337 | 0.03482339 | 4.9e-07 |
| compat_eps0.01_t10.0 | 0.01885618 | 0.01885624 | 3.0e-06 |
| compat_eps0.1_t0.5 | 0.35680774 | 0.35680774 | 0.0e+00 |
| compat_eps0.1_t1.0 | 0.33063745 | 0.33063745 | 1.7e-16 |

## 7. (e) Code audit: bugs and methodological issues (no code was changed)

Severity: **A** = can change a reported conclusion; **B** = bias/ambiguity that should be disclosed or ablated; **C** = minor. Line numbers refer to the files in `D:\Videos\KazNPU\Symmetry_PINN_2026\code`.

### A. High impact

* **A1. FD reference: first-order Neumann rows** - `common.py:120-121` (`ab[1,0]=1+lam; ab[1,-1]=1+lam` with off-diagonals `-lam`) discretise `u_x=0` as `u_1-u_0=0` (first order; equivalent to dropping half the diffusive flux at the wall node). Consistent ghost-node rows are `1+2lam`, `-2lam`. Errors vs trusted: Section 4. This reference feeds every table (all PINN errors, the Chebyshev baseline, the eps/K/T sweeps) and the inverse-problem data (`inverse_problem.py:25`). Backward Euler (`common.py:114-122`) is additionally only first order in dt (negligible against the wall-row error).
* **A2. Chebyshev baseline scored against that reference with a mismatched norm** - `chebyshev_spectral.py:84-86`, `run_experiments.py:86-87`, `stress_eps_sweep.py:80-82`: FD data are `np.interp`-ed onto the (non-uniform, N+1) Chebyshev nodes and `l2_error(u, u_ref, dx_FD)` sums `e_j^2*dx_FD` over those nodes (total "length" 81*0.004=0.32 instead of 2; 401*0.003=1.2 in the sweep). The stored numbers are dominated by the reference's own error (Section 4.3): "no improvement with N, dt", the Linf floor 5.5e-4, and the eps-sweep decrease with eps are artifacts. The corrected Chebyshev errors are 1.7e-4 / 8.4e-5 / 4.2e-5 (N=40/60/80). The baseline is only first order in time (backward Euler, `chebyshev_spectral.py:32-37`), so "classical collocation is 0.16 s" is not its best accuracy per cost: dt=5e-5 gives 4.2e-6 in 0.46 s (my LU replica), a certified 1e-10 reference needs 0.2 s.
* **A3. Boundary-condition loss effectively switched off, and it conflicts with the IC fit** - `spectral_pinn.py:72-74` divides the boundary derivative vectors by `K^2`, and `:160-161` squares the residual: the loss weight on `u_x(+-1)^2` is `lam_bc/K^4` = 9.5e-7 (K=32), 6e-8 (K=64), 2.3e-10 (K=256). The IC (`:80`, `c0` = LS fit of `u0`) has `u_x(+-1)=-+2`; the resulting BC penalty at `c0` is 7.6e-6 at K=32, 4.8e-7 at K=64, 3.0e-8 at K=128 (table below), against a final total training loss of 4.1e-5 at K=32 (`spectral_pinn_energy.json`, `loss_history`). The trained K=32 network (seed 42) has wall slopes -0.38, +0.24, +0.13, -0.26 at t=0.2, 0.5, 0.8, 1.0 (true: 0). Wall errors are not what limits it (|x|>0.95: Linf 3.6e-3 at t=0.2, 3.5e-4 at t=1.0; interior Linf 4.4e-3 / 3.8e-3 at |x|~0.67-0.76), but the Neumann condition is not actually enforced, the incompatible IC is fitted exactly at t=0 (`:157-158`), and for K>=64 the BC term is ~inert, which confounds the K-sweep. The vanilla PINN has a separate treatment (`vanilla_pinn.py:53-57,77,79`: BC weight 10, applied also at t=0 where it contradicts the IC at (+-1,0)).

| K | u_x(+1) of IC fit c0 | L_bc(c0) as coded (scaled by 1/K^4) | scale factor 1/K^4 |
|---|---|---|---|
| 16 | -2 | 1.22e-04 | 1.53e-05 |
| 32 | -2 | 7.63e-06 | 9.54e-07 |
| 64 | -2 | 4.77e-07 | 5.96e-08 |
| 128 | -2.001 | 2.98e-08 | 3.73e-09 |
| 256 | -1.817e+04 | 1.54e-01 | 2.33e-10 |

| t | u_x(+1), Spectral-PINN seed 42 | u_x(-1) | u_x(+1) of K=32 projection of trusted solution (true value 0) |
|---|---|---|---|
| 0.2 | -0.376 | 0.381 | -0.056 |
| 0.5 | 0.241 | -0.222 | 0.069 |
| 0.8 | 0.132 | -0.106 | 0.028 |
| 1.0 | -0.261 | 0.231 | 0.018 |

* **A4. Ill-conditioned equispaced grid for IC fit and PDE collocation at large K** - `spectral_pinn.py:60,65,77,80`: `n_grid=300` equispaced points for every K; `fit_chebyshev_coeffs` (`common.py:92-96`) is a least-squares fit of a degree-K Chebyshev series on them. Rule of thumb: stable polynomial fitting on n equispaced points needs degree of order sqrt(n) (~2 sqrt(n) = 35 for n=300). Measured:

| K | cond(Vandermonde 300 equispaced pts) | IC-fit error at the 300 pts | IC-fit sup error (40001 pts) | wall slope of IC fit (true -2) |
|---|---|---|---|---|
| 16 | 3.8e+00 | 6.2e-11 | 6.2e-11 | -2 |
| 32 | 4.0e+00 | 2.2e-15 | 2.9e-15 | -2 |
| 48 | 6.5e+00 | 2.2e-15 | 2.4e-15 | -2 |
| 64 | 4.1e+01 | 1.4e-15 | 4.9e-15 | -2 |
| 96 | 8.9e+04 | 3.0e-15 | 8.7e-12 | -2 |
| 128 | 1.4e+10 | 2.2e-15 | 8.7e-07 | -2.001 |
| 192 | 3.8e+15 | 7.1e-15 | 1.1e+00 | -5225 |
| 256 | 7.5e+15 | 1.4e-14 | 1.2e+00 | -1.817e+04 |

  For K>=192 the "IC fit" `c0` matches `u0` at the 300 points (1e-14) but is O(1) wrong (1.1-1.2) between them, with wall slopes -5.2e3 / -1.8e4; the warm-up (`:139-146`) then trains c(t) to this c0. The K=256 results in `stress_k_sweep.py` (L2 1.51, 1.47) are therefore confounded by (and plausibly explained by) the initialisation/collocation conditioning rather than demonstrably by "optimization difficulty", and the K=128 (cond 1e10) and K=64 rows are likewise confounded (the PDE residual is only sampled at 300 equispaced points, which cannot control a degree-128/256 polynomial between them; for the T-sweep this is visible directly in the K-floor table, 4.5). The abstract/Section "stress test" claim that accuracy degrades with K "exposing an optimization bottleneck absent in classical spectral methods" is not supported by this experiment. (K=64 at 1x vs 2x budget: L2 6.2e-2 -> 2.3e-2 (`data/stress_k_sweep.json`), so even there the budget matters; no 4x/8x data for K=64.)
* **A5. The energy penalty is never active, and the violation test is vacuous** - `spectral_pinn.py:165-171`: `lam_energy*sum(clamp(dE,0)^2)` over 11 checkpoints (`:133`; `dE`= difference of the network's own energies). The warm-up (`:139-146`) starts from `c(t)=c0` for all t, hence `dE=0` initially; afterwards the true dynamics decrease E, so the clamp is zero. `data/spectral_pinn_energy.json: energy_history` (every 50th epoch, 120 samples) is nonzero only at the first sample (7.6e-8) and exactly 0 afterwards; the runs with/without penalty agree to <1e-3 relative (L2 1.4517e-3 vs 1.4505e-3 at t=1; seed means 1.3777e-3 vs 1.3785e-3) and their E curves to 4.5e-7. The "zero energy violations" count (`run_experiments.py:65`, `stress_time_sweep.py:57`: increase > 1e-10 between consecutive evaluation points, 51 points of spacing 0.02 at T=1; 20T+1 points of spacing 0.05 for T>1) cannot trigger in the T=1 regime (Section 6). The penalty weight actually used is `lam_energy=50` (`spectral_pinn.py:125`), 10 checkpoints, never varied; there is no ablation in the repo.
* **A6. The benchmark does not stress eps (or diffusion)** - Section 5. The eps-sweep (`stress_eps_sweep.py:26`) and the claims about interface sharpness at small eps are not tested by a T=1 problem whose solution is within 3e-3 of the eps=0 ODE solution. The inverse problem at eps=0.01 is nearly unidentifiable (5.5).

### B. Bias / disclosure

* **B1. Inverse data generated from the faulty reference** - `inverse_problem.py:25-31` (`reference_solution_fd` + `np.interp`): RMS data error 3.3e-4 (Linf 2.4e-3) at eps=0.01, 6.2e-4 (Linf 3.4e-3) at eps=0.1 (5.5).
* **B2. Inconsistent / untuned loss weights between the main inverse run and the sweep; no convergence record** - main: `inverse_problem.py:36-37` `lam_pde=10, lam_data=5, lam_ic=10, lam_bc=1 (K^-4-scaled), eps_init=0.03` (3x the truth); sweep: `server_sweep_inverse.py:25-29` `lam_pde=50, lam_data=2, eps_init=2*eps_true`; both give ~+38-40%. The data term averages 200 values (5 times x 40 points) while the PDE term averages 100x300 residuals (`:55`, `spectral_pinn.py:65`); the loss change caused by a 40% change of eps is `lam_data * RMS^2 ~ 5*(7e-4)^2 = 2.5e-6` at eps=0.01. The scalar `log_eps` gets 10x lr (`:49-52`) but the same cosine schedule (`:53`) which decays to 0, so eps freezes in the last part of training; `run_experiments.py:100` stores `eps_history[::50]` but no convergence criterion is used. One noise realization only (`make_observations(..., seed=1)`, `:39`). The main-text inverse results are not archived: `run_experiments.py:102` writes `inverse_problem.json`, which does not exist in `data/` (only `inverse_bias_sweep.json`).
* **B3. What the bias is (and is not)** - (no training; `code_checks.py`) a least-squares estimator of eps^2 from the PDE residual, evaluated with the K-projected trusted solution on the repo's own collocation grid (100 t-points from 1e-4, 300 x-points), returns eps_hat/eps_true as follows (no systematic +40%; scatter 0.69-1.39 from the unresolved wall layer at small eps and large K):

| eps_true | K=16 | K=32 | K=64 | K=96 |
|---|---|---|---|---|
| 0.005 | 1.093 | 1.388 | 0.956 | 1.167 |
| 0.01 | 1.109 | 1.091 | 0.910 | 0.686 |
| 0.02 | 1.076 | 1.004 | 0.964 | 0.863 |
| 0.05 | 0.989 | 1.020 | 1.174 | 1.246 |
| 0.1 | 0.999 | 1.055 | 1.142 | 1.103 |

  The sweep results (`data/inverse_bias_sweep.json`) give eps_hat/eps_true = 1.28 +- 0.24 (K=16), 1.400 +- 0.034 (K=32), 1.372 +- 0.014 (K=64), 1.409 +- 0.021 (K=96); pooled K>=32 (n=15): **1.394 +- 0.028, vs sqrt(2) = 1.414; mean of squares 1.944 ~ 2**, for eps_true from 0.005 to 0.1 and two different initialisations (3x in the main run, 2x in the sweep). A constant factor ~sqrt(2) in eps (2 in eps^2) over a 20x range of eps and K from 16 to 96 (K>=32 pooled) is a systematic effect of the joint optimisation (weights/IC-BC conflict/data design), not a truncation artifact and not noise. I did not run ablations; this is evidence for where not to look, not a diagnosis.
* **B4. Soft IC; no input scaling at long T** - `spectral_pinn.py:157-158`: IC enforced by a coefficient-space MSE with weight 10 at one time point. At longer horizons the network drifts: E(0) of the trained network is 0.400457 (T=1), 0.401214 (T=5), 0.403848 (T=10) vs exact 0.400436 (`stress_time_sweep.json`: `E_first`). `t` enters the tanh MLP unscaled (`stress_time_sweep.py:23`, `t in [0,10]`). Warm-up (`:139-146`) fixes the start at c(t)=c0.
* **B5. Collocation in time** - `spectral_pinn.py:131`: fixed 100 equispaced points `linspace(1e-4,T,100)`, never resampled (`stress_time_sweep.py:29`: 100*T points, so the density is constant); the interval (0, 1e-4) and the thin early wall layer (eps sqrt(t) < K=32 spacing 4.8e-3 for t<0.23) are not represented; energy checkpoints are 11 points (`:133`). (Empirically the wall region is not what limits the K=32 error, A3.)
* **B6. Multi-seed spread is not a robustness measure** - `server_sweep_forward.py:18`: seeds change only the network initialisation; the warm-up (`spectral_pinn.py:139-146`) maps every seed to the same constant-in-time c0 and collocation is deterministic. Ratio statistics: vanilla L2 0.0674 +- 0.0055, Spectral 0.001379 +- 7.7e-5 (n=5); **ratio of means 48.9** (bootstrap 95% CI 45.0-52.6, 20000 resamples) = the abstract's "~49"; seed-42 single-run ratio 46.5 (0.067544/0.001451) = the Conclusions' "~47"; per-seed ratios 46.5, 50.2, 54.7, 50.9, 42.5. Wall-clock ratio 919.3/169.4 = 5.43 (abstract "5.4x"), measured on the parallel server with 2 threads/process; the Spectral-PINN time is inflated by the per-column autograd loop (`spectral_pinn.py:97-100`: K+1 backward passes) and the vanilla one by second derivatives on 6400 points - an implementation comparison, not a complexity statement.
* **B7. Vanilla baseline strength** - `vanilla_pinn.py:35-37,43-47,79`: one configuration (3x64 tanh MLP, Adam 8000 epochs, lr 1e-3 cosine to 0, no input scaling, no L-BFGS, no causal/adaptive weights, lam_ic=lam_bc=10, uniform 80x80 collocation); relative L2 error 4.6-9.2% (t=0.2-1.0), i.e. 25-50x worse than the closed-form ODE solution (5.4). The "factor ~49" is relative to this untuned baseline at equal epochs, not equal budget/tuning.
* **B8. 2D symmetry experiment** - `experiment_2d_symmetry.py:81-85,166-178`: the cos-cos basis is even by construction, so S(t)=0 to rounding (1e-31) is a property of the basis, not of training; the vanilla comparator is given `sin x, sin y` inputs (`:203`) so it *can* break the symmetry - a vanilla PINN with cos(x), cos(y) inputs would also be exactly symmetric. There is no accuracy metric vs a reference solution in 2D at all; learning rates differ (`:113` 2e-3 vs `:211` 1e-3), collocation sizes differ (48^2 grid x 40 times vs 24^2 x 15, `:112,210`); `EPS2D=0.3`, amplitude 0.4 (`:36,46`) is a smooth decaying regime without interfaces.
* **B9. Long-time stress test** - `stress_time_sweep.py:23-57`: besides the K=32 representation floor (4.5: L2 >= 0.135-0.157 at t>=5, max|fit| up to 1.18 > 1, so a K=32 series cannot respect the bound principle there), the network's E_last (0.249 at T=5, 0.263 at T=10) is 4.8x the true E(5)=0.0514 and 14x the true plateau value 0.0189 (T=10): a monotone decrease of a wrong energy.

### C. Minor

* **C1.** `common.py:130-131` `l2_error = sqrt(sum e^2 dx)`: endpoints carry full weight; equals the continuous L2 to <1% for smooth errors, overstates it up to 2.8x for wall-layer errors (4.1).
* **C2.** Linear `np.interp` of the FD reference everywhere (`run_experiments.py:33,53,86`; `vanilla_pinn.py:105`; `server_sweep_forward.py:31,48`; `stress_eps_sweep.py:40,60,80`; `stress_k_sweep.py:37`; `stress_time_sweep.py:47`; `inverse_problem.py:29`): floor L2_200 1.3e-5 .. 3.3e-5, Linf up to 1.2e-4 (Nx=500). Irrelevant for PINN errors of 1e-3, comparable to the corrected Chebyshev errors.
* **C3.** The reference differs between experiments: `reference_nx_for_eps` (`stress_common.py:7-10`) gives Nx=668 at eps=0.01 in the eps-sweep vs Nx=500 elsewhere (FD-ref error 6.8e-5 / 4.3e-4 vs 8.9e-5 / 5.6e-4), so the "same" eps=0.01 Chebyshev row differs between Table 1 (8.5e-5) and the sweep (1.17e-4).
* **C4.** `stress_k_sweep.py:63`: `baseline_K32` (1.38e-3, 3.32e-3) is typed in, not computed (equals the 5-seed means 1.3785e-3, 3.3193e-3).
* **C5.** `inverse_problem.py:45`: `SpectralPINN(eps_true, ...)` stores the true eps in the model object (unused in the loss; harmless, but leak-prone).
* **C6.** Maximum bound principle: the trusted solutions satisfy |u|<=1; stored Spectral-PINN predictions lie in [-0.9971, 0.1493] (t<=1); the vanilla PINN at t=1 reaches -1.00079 (|u|>1 by 8e-4); at long times the K=32 series cannot satisfy |u|<=1 (max|fit| 1.17-1.18 at t>=5). The experiments never check it.
* **C7.** Learning rates / schedules (no tuning or ablation reported): Spectral-PINN Adam 1e-3 (`spectral_pinn.py:123`), warm-up Adam 1e-3 (`:139`), cosine to 0; inverse Adam 2e-3 for the network and 2e-2 for `log_eps` (`inverse_problem.py:36,49-52`); 2D Fourier-PINN 2e-3 (`experiment_2d_symmetry.py:113`) vs 2D vanilla 1e-3 (`:211`); vanilla 1e-3 (`vanilla_pinn.py:37`). Epoch counts: Spectral 6000 (+1000 warm-up), vanilla 8000, inverse 4000.

## 8. The module `refsol.py`

```
import refsol
refsol.u0_orig(x), refsol.u0_compat(x)            # x^2 cos(pi x);  (cos(pi x) - cos^2(pi x))/2
res = refsol.reference(eps, ic, T, times, x_eval)  # ic: 'orig' | 'compat' | callable; dict time -> values on x_eval
res, err_est = refsol.reference_certified(...)     # additionally returns the N-vs-1.5N difference
refsol.l2_repo(e, x), refsol.l2_trapz(e, x), refsol.linf(e), refsol.bary_interp(xn, un, x), refsol.energy_cc(eps, xn, un)
refsol.self_test()                                 # or: python refsol.py   (6-22 s)
```

Self-test (all PASS, 6.0 s on an idle machine, 22 s under load): IC derivative values (u0'(-1)=2, u0'(1)=-2, u0c'(+-1)=0, compat identities); differentiation matrix, barycentric interpolation, Clenshaw-Curtis weights; exact solution of the linearised problem (`a cos(k pi x) e^{(1-eps^2 k^2 pi^2) t}`, relative error 6.5e-12 / 3.0e-12 for k=1,2); N-convergence for both ICs; independent FD-Richardson comparison; certified-wrapper estimate for eps=0.02, T=10 (1.4e-11) and the long-time energy equal to the two-interface value `2(2 sqrt2/3) eps` to < 1e-4; half- vs full-domain; BC, bound principle, parity, energy monotonicity; Radau vs BDF. Typical cost: 0.2 s per `reference(...)` call with 3 output times at T=1 (matrices cached), 10 s for T=10, eps=0.01 (N=1024).
