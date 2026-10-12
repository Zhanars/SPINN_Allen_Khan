"""
Spectral-PINN hybrid for the 1D Allen-Cahn equation.

Idea (grant task 3 / "Spectral-PINN гибридті архитектурасы"): the neural
network does not output u(x,t) directly. Instead it outputs a
time-dependent vector of Chebyshev coefficients  c(t) = (c_0(t),...,c_K(t)),
and the field is reconstructed as

    u(x,t) = sum_{k=0}^{K} c_k(t) T_k(x).

Spatial derivatives are obtained *exactly* by a fixed linear operator in
coefficient space (Dcoef, see common.py) rather than by automatic
differentiation, which removes the spectral bias / sharp-interface
difficulty that plagues vanilla space-time PINNs at small eps. Automatic
differentiation is only used in the time direction (dc/dt).

An optional "energy-constrained" loss term implements the discrete energy
dissipation theorem (grant task 3): the Ginzburg-Landau free energy
evaluated on the network's own trajectory is penalized whenever it
increases between two checkpoints.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import time
import numpy as np
import torch
import torch.nn as nn
from common import (u0_fn, f_ac, F_ac, T_FINAL,
                     cheb_coeff_derivative_matrix, chebyshev_vandermonde,
                     cheb_boundary_deriv_vectors, fit_chebyshev_coeffs)

torch.set_default_dtype(torch.float64)


class CoeffNet(nn.Module):
    """t -> Chebyshev coefficient vector c(t) in R^{K+1} (or K-1 free coefficients with the hard Neumann basis)."""
    def __init__(self, K, hidden=(64, 64), n_out=None):
        super().__init__()
        dims = [1, *hidden, K + 1 if n_out is None else n_out]
        layers = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            if i < len(dims) - 2:
                layers.append(nn.Tanh())
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, t):
        return self.net(t)


def grad(y, x):
    return torch.autograd.grad(y, x, grad_outputs=torch.ones_like(y), create_graph=True)[0]


class SpectralPINN:
    def __init__(self, eps, K=32, hidden=(64, 64), n_grid=300, seed=42, u0=None, grid="uniform", hard_bc=False):
        """grid='uniform': 300 equispaced points (as submitted); grid='cgl': Chebyshev-Gauss-Lobatto nodes
        (N = max(n_grid, 4K+8)) with Clenshaw-Curtis weights.  hard_bc=True: the network outputs K-1 free
        coefficients b and c = M b spans exactly the polynomials with u_x(+-1)=0 (basis T_k - k^2/(k+2)^2 T_{k+2})."""
        torch.manual_seed(seed)
        self.eps = eps
        self.K = K
        self.grid = grid
        self.hard_bc = hard_bc
        if hard_bc:
            M = np.zeros((K + 1, K - 1))
            for k in range(K - 1):
                M[k, k] = 1.0
                M[k + 2, k] = -(k ** 2) / (k + 2) ** 2
            self.M = torch.tensor(M)
        self.net = CoeffNet(K, hidden, n_out=(K - 1 if hard_bc else None))
        if grid == "cgl":
            N = max(n_grid, 4 * K + 8)
            self.x_grid = np.cos(np.pi * np.arange(N + 1) / N)[::-1].copy()
            from chebyshev_spectral import clenshaw_curtis_weights
            self.w_np = clenshaw_curtis_weights(N)[::-1].copy()
            self.dx = None
        else:
            self.x_grid = np.linspace(-1, 1, n_grid)
            self.dx = self.x_grid[1] - self.x_grid[0]
            self.w_np = None
        Dc_np = cheb_coeff_derivative_matrix(K)
        Tm_np = chebyshev_vandermonde(self.x_grid, K)
        vp_np, vm_np = cheb_boundary_deriv_vectors(K)
        # Normalize boundary vectors (entries grow like k^2) so the BC loss
        # is O(1)-scaled and does not dominate the PDE/IC gradients.
        bc_scale = max(1.0, K**2)
        vp_np = vp_np / bc_scale
        vm_np = vm_np / bc_scale
        self.Dc = torch.tensor(Dc_np)
        self.Dc2 = torch.tensor(Dc_np @ Dc_np)
        self.Tm = torch.tensor(Tm_np)             # (n_grid, K+1)
        self.v_plus = torch.tensor(vp_np)
        self.v_minus = torch.tensor(vm_np)
        self.u0 = u0_fn if u0 is None else u0
        if grid == "cgl":
            sw = np.sqrt(self.w_np)
            c0, *_ = np.linalg.lstsq(Tm_np * sw[:, None], self.u0(self.x_grid) * sw, rcond=None)
            self.wq = torch.tensor(self.w_np)
        else:
            c0 = fit_chebyshev_coeffs(self.u0(self.x_grid), self.x_grid, K)
        self.c0 = torch.tensor(c0)

    def coeffs(self, t):
        out = self.net(t)
        return out @ self.M.T if self.hard_bc else out

    def pde_loss(self, res):
        if self.grid == "cgl":
            return ((res ** 2) * self.wq).sum(dim=-1).mean() / 2.0
        return (res ** 2).mean()

    def field(self, c):
        """c: (batch, K+1) -> u: (batch, n_grid)"""
        return c @ self.Tm.T

    def pde_residual(self, t_col):
        t_col = t_col.clone().requires_grad_(True)
        c = self.coeffs(t_col)                      # (B, K+1)
        # dc_k/dt for every coefficient k (independent columns), obtained
        # one column at a time since autograd.grad with grad_outputs=1
        # only returns the *sum* over output dimensions of the VJP.
        cols = []
        for k in range(self.K + 1):
            gk = torch.autograd.grad(c[:, k], t_col, grad_outputs=torch.ones_like(c[:, k]),
                                      create_graph=True, retain_graph=True)[0]
            cols.append(gk)
        dcdt = torch.cat(cols, dim=1)                # (B, K+1)
        u = self.field(c)                            # (B, n_grid)
        u_t = self.field(dcdt)
        u_xx_coeff = c @ self.Dc2.T
        u_xx = self.field(u_xx_coeff)
        res = u_t - self.eps**2 * u_xx + f_ac(u)
        return res, c

    def bc_residual(self, c):
        r_plus = c @ self.v_plus
        r_minus = c @ self.v_minus
        return r_plus, r_minus

    def energy(self, c):
        """Physical Ginzburg-Landau energy via trapezoid rule on x_grid."""
        u = self.field(c)
        ux_coeff = c @ self.Dc.T
        ux = self.field(ux_coeff)
        integrand = self.eps**2 / 2 * ux**2 + F_ac(u)
        if self.grid == "cgl":
            return (integrand * self.wq).sum(dim=-1)
        return torch.trapezoid(integrand, dx=self.dx, dim=-1)


def train_spectral_pinn(eps, K=32, epochs=6000, nt_col=100, lr=1e-3,
                         lam_ic=10.0, lam_bc=1.0, energy_constrained=False,
                         lam_energy=50.0, n_checkpoints=10, seed=42,
                         warmup_epochs=1000, verbose=True, u0=None, T=T_FINAL, grid="uniform", hard_bc=False, lbfgs_iters=0):
    model = SpectralPINN(eps, K=K, seed=seed, u0=u0, grid=grid, hard_bc=hard_bc)
    opt = torch.optim.Adam(model.net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    t_col = torch.linspace(1e-4, T, nt_col).reshape(-1, 1)
    t_ic = torch.zeros(1, 1)
    t_ckpt = torch.linspace(0, T, n_checkpoints + 1).reshape(-1, 1)

    # Warm-up: pretrain the network to output the (constant-in-time) IC
    # coefficients everywhere. This gives the PDE-constrained optimization
    # a sane starting point instead of an arbitrary random trajectory,
    # analogous to the per-step warm-up used in discrete SAV-PINN schemes.
    warmup_opt = torch.optim.Adam(model.net.parameters(), lr=1e-3)
    t_warmup = torch.linspace(0, T, nt_col).reshape(-1, 1)
    for _ in range(warmup_epochs):
        warmup_opt.zero_grad()
        c_pred = model.coeffs(t_warmup)
        loss_w = ((c_pred - model.c0)**2).mean()
        loss_w.backward()
        warmup_opt.step()
    if verbose:
        print(f"  [warmup] final loss = {loss_w.item():.3e}")

    t0 = time.time()
    history = {"loss": [], "l_pde": [], "l_ic": [], "l_bc": [], "l_energy": []}
    for ep in range(epochs):
        opt.zero_grad()
        res, c_col = model.pde_residual(t_col)
        L_pde = model.pde_loss(res)

        c_ic = model.coeffs(t_ic)
        L_ic = ((c_ic - model.c0)**2).mean()

        r_plus, r_minus = model.bc_residual(c_col)
        L_bc = (r_plus**2).mean() + (r_minus**2).mean()

        loss = L_pde + lam_ic * L_ic + lam_bc * L_bc

        L_energy = torch.tensor(0.0)
        if energy_constrained:
            c_ck = model.coeffs(t_ckpt)
            E_ck = model.energy(c_ck)
            dE = E_ck[1:] - E_ck[:-1]
            L_energy = (torch.clamp(dE, min=0.0)**2).sum()
            loss = loss + lam_energy * L_energy

        loss.backward()
        opt.step()
        sched.step()
        history["loss"].append(loss.item())
        history["l_pde"].append(L_pde.item())
        history["l_ic"].append(L_ic.item())
        history["l_bc"].append(L_bc.item())
        history["l_energy"].append(float(L_energy.detach()) if torch.is_tensor(L_energy) else float(L_energy))
        if verbose and (ep % 1000 == 0 or ep == epochs - 1):
            print(f"  ep {ep:>5} loss={loss.item():.3e} pde={L_pde.item():.3e} "
                  f"ic={L_ic.item():.3e} bc={L_bc.item():.3e} energy={float(L_energy):.3e}")
    if lbfgs_iters > 0:
        # full-batch L-BFGS refinement of the same loss (the collocation grid is fixed, so the loss is deterministic)
        optL = torch.optim.LBFGS(model.net.parameters(), lr=1.0, max_iter=lbfgs_iters, max_eval=int(1.25 * lbfgs_iters),
                                 history_size=50, line_search_fn="strong_wolfe", tolerance_grad=1e-14, tolerance_change=1e-16)

        def closure():
            optL.zero_grad()
            res_, c_col_ = model.pde_residual(t_col)
            l = model.pde_loss(res_) + lam_ic * ((model.coeffs(t_ic) - model.c0)**2).mean()
            rp_, rm_ = model.bc_residual(c_col_)
            l = l + lam_bc * ((rp_**2).mean() + (rm_**2).mean())
            if energy_constrained:
                E_ = model.energy(model.coeffs(t_ckpt))
                l = l + lam_energy * (torch.clamp(E_[1:] - E_[:-1], min=0.0)**2).sum()
            l.backward()
            return l
        optL.step(closure)
        history["lbfgs_final_loss"] = float(closure())
    elapsed = time.time() - t0
    return model, history, elapsed


def eval_on_grid(model, t_val, x_grid=None):
    with torch.no_grad():
        c = model.coeffs(torch.tensor([[t_val]]))
        if x_grid is None:
            return model.field(c).numpy().flatten()
        Tm = torch.tensor(chebyshev_vandermonde(x_grid, model.K))
        return (c @ Tm.T).numpy().flatten()


if __name__ == "__main__":
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]

    print("=== Spectral-PINN (no energy constraint) ===")
    model, hist, elapsed = train_spectral_pinn(eps, K=32, epochs=6000, energy_constrained=False)
    u_pred = eval_on_grid(model, 1.0, x_ref)
    print(f"L2={l2_error(u_pred, ref_snaps[1.0], dx_ref):.4e} "
          f"Linf={linf_error(u_pred, ref_snaps[1.0]):.4e} time={elapsed:.0f}s")

    print("=== Spectral-PINN + energy-constrained loss ===")
    model_e, hist_e, elapsed_e = train_spectral_pinn(eps, K=32, epochs=6000, energy_constrained=True)
    u_pred_e = eval_on_grid(model_e, 1.0, x_ref)
    print(f"L2={l2_error(u_pred_e, ref_snaps[1.0], dx_ref):.4e} "
          f"Linf={linf_error(u_pred_e, ref_snaps[1.0]):.4e} time={elapsed_e:.0f}s")
