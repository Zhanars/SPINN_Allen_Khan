"""Baseline: standard (vanilla) space-time PINN for the 1D Allen-Cahn equation."""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import time
import numpy as np
import torch
import torch.nn as nn
from common import u0_fn, f_ac, T_FINAL

torch.set_default_dtype(torch.float64)


class MLP(nn.Module):
    def __init__(self, dims):
        super().__init__()
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

    def forward(self, xt):
        return self.net(xt)


def grad(y, x):
    return torch.autograd.grad(y, x, grad_outputs=torch.ones_like(y), create_graph=True)[0]


def train_vanilla_pinn(eps, seed=42, epochs=8000, nx_col=80, nt_col=80,
                        n_ic=200, n_bc=200, lr=1e-3, lam_ic=10.0, lam_bc=10.0,
                        arch=(2, 64, 64, 64, 1), verbose=True, u0=None, T=T_FINAL, lbfgs_iters=0):
    torch.manual_seed(seed)
    net = MLP(list(arch))
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    x_col = torch.linspace(-1, 1, nx_col)
    t_col = torch.linspace(0, T, nt_col)
    Xc, Tc = torch.meshgrid(x_col, t_col, indexing="ij")
    Xc = Xc.reshape(-1, 1).requires_grad_(True)
    Tc = Tc.reshape(-1, 1).requires_grad_(True)

    x_ic = torch.linspace(-1, 1, n_ic).reshape(-1, 1)
    t_ic = torch.zeros_like(x_ic)
    u_ic = torch.tensor((u0_fn if u0 is None else u0)(x_ic.numpy()))

    t_bc = torch.linspace(0, T, n_bc).reshape(-1, 1)
    xl_bc = -torch.ones_like(t_bc).requires_grad_(True)
    xr_bc = torch.ones_like(t_bc).requires_grad_(True)
    t_bc_l = t_bc.clone().requires_grad_(True)
    t_bc_r = t_bc.clone().requires_grad_(True)

    t0 = time.time()
    history = []
    for ep in range(epochs):
        opt.zero_grad()
        u = net(torch.cat([Xc, Tc], dim=1))
        u_x = grad(u, Xc)
        u_xx = grad(u_x, Xc)
        u_t = grad(u, Tc)
        res = u_t - eps**2 * u_xx + f_ac(u)
        L_pde = (res**2).mean()

        u_ic_pred = net(torch.cat([x_ic, t_ic], dim=1))
        L_ic = ((u_ic_pred - u_ic)**2).mean()

        u_l = net(torch.cat([xl_bc, t_bc_l], dim=1))
        u_r = net(torch.cat([xr_bc, t_bc_r], dim=1))
        ux_l = grad(u_l, xl_bc)
        ux_r = grad(u_r, xr_bc)
        L_bc = (ux_l**2).mean() + (ux_r**2).mean()

        loss = L_pde + lam_ic * L_ic + lam_bc * L_bc
        loss.backward()
        opt.step()
        sched.step()
        if verbose and (ep % 1000 == 0 or ep == epochs - 1):
            print(f"  ep {ep:>5} loss={loss.item():.3e} pde={L_pde.item():.3e} "
                  f"ic={L_ic.item():.3e} bc={L_bc.item():.3e}")
        history.append(loss.item())
    if lbfgs_iters > 0:
        optL = torch.optim.LBFGS(net.parameters(), lr=1.0, max_iter=lbfgs_iters, max_eval=int(1.25 * lbfgs_iters),
                                 history_size=50, line_search_fn="strong_wolfe", tolerance_grad=1e-14, tolerance_change=1e-16)

        def closure():
            optL.zero_grad()
            u = net(torch.cat([Xc, Tc], dim=1))
            u_x = grad(u, Xc); u_xx = grad(u_x, Xc); u_t = grad(u, Tc)
            L_pde = ((u_t - eps**2 * u_xx + f_ac(u))**2).mean()
            L_ic = ((net(torch.cat([x_ic, t_ic], dim=1)) - u_ic)**2).mean()
            u_l = net(torch.cat([xl_bc, t_bc_l], dim=1)); u_r = net(torch.cat([xr_bc, t_bc_r], dim=1))
            L_bc = (grad(u_l, xl_bc)**2).mean() + (grad(u_r, xr_bc)**2).mean()
            l = L_pde + lam_ic * L_ic + lam_bc * L_bc
            l.backward()
            return l
        optL.step(closure)
    elapsed = time.time() - t0
    return net, history, elapsed


def eval_on_grid(net, x_grid, t_val):
    with torch.no_grad():
        xt = torch.stack([torch.tensor(x_grid), torch.full_like(torch.tensor(x_grid), t_val)], dim=1)
        return net(xt).numpy().flatten()


if __name__ == "__main__":
    from common import reference_solution_fd, l2_error, linf_error
    eps = 0.01
    x_ref, ref_snaps = reference_solution_fd(eps)
    dx_ref = x_ref[1] - x_ref[0]
    net, hist, elapsed = train_vanilla_pinn(eps, epochs=8000)
    x_test = np.linspace(-1, 1, 200)
    u_pred = eval_on_grid(net, x_test, 1.0)
    u_ref_i = np.interp(x_test, x_ref, ref_snaps[1.0])
    print(f"Vanilla PINN: L2={l2_error(u_pred,u_ref_i,x_test[1]-x_test[0]):.4e} "
          f"Linf={linf_error(u_pred,u_ref_i):.4e} time={elapsed:.0f}s")
