"""
2D Allen-Cahn on the periodic square [-pi,pi]^2 (grant task 4: "1D and 2D
numerical experiments"), used to make a *quantitative* symmetry argument
appropriate for the journal Symmetry.

Initial condition u0(x,y) = 0.4*cos(2x)*cos(2y) is invariant under the
dihedral symmetries of the square: x -> -x, y -> -y, and (x,y) -> (y,x).
Because the Allen-Cahn operator (isotropic Laplacian + the odd
double-well nonlinearity f(u)=u^3-u) commutes with these reflections, the
exact solution stays symmetric for all t>0.

- Fourier-Spectral-PINN: represents u(x,y,t) = sum_{m,n=0}^{M}
  c_{mn}(t) cos(mx)cos(ny). Because every basis function is already even
  in x and in y, u(-x,y,t)=u(x,y,t) and u(x,-y,t)=u(x,y,t) HOLD EXACTLY,
  to machine precision, regardless of training quality. The Laplacian is
  diagonal in mode space: Delta[cos(mx)cos(ny)] = -(m^2+n^2)cos(mx)cos(ny).

- Vanilla PINN: takes (cos x, sin x, cos y, sin y, t) as input (so
  periodicity is automatic) and predicts u directly via an MLP; nothing
  in its architecture forces the reflection symmetries, so any violation
  is purely a training/approximation artifact.

We measure the reflection-symmetry violation
    S(t) = mean_{x,y}( (u(x,y,t)-u(-x,y,t))^2 + (u(x,y,t)-u(x,-y,t))^2 )
for both methods.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import time
import numpy as np
import torch
import torch.nn as nn

torch.set_default_dtype(torch.float64)

EPS2D = 0.3
T2D = 1.0
M_MODES = 8


def f_ac(u):
    return u**3 - u


def u0_2d(x, y):
    return 0.4 * np.cos(2 * x) * np.cos(2 * y)


# ------------------------------------------------------------------
# Fourier-Spectral-PINN (cos-cos basis, exact reflection symmetry)
# ------------------------------------------------------------------
class FourierCoeffNet(nn.Module):
    def __init__(self, M, hidden=(64, 64)):
        super().__init__()
        n_out = (M + 1) * (M + 1)
        dims = [1, *hidden, n_out]
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
        self.M = M

    def forward(self, t):
        out = self.net(t)
        return out.reshape(-1, self.M + 1, self.M + 1)


class FourierSpectralPINN2D:
    def __init__(self, eps, M=M_MODES, n_grid=48, seed=42):
        torch.manual_seed(seed)
        self.eps = eps
        self.M = M
        self.net = FourierCoeffNet(M)
        m = np.arange(M + 1)
        self.lap_eig = torch.tensor(-(m[:, None]**2 + m[None, :]**2), dtype=torch.float64)
        grid = np.linspace(-np.pi, np.pi, n_grid, endpoint=False)
        self.x_grid = grid
        self.y_grid = grid
        self.Cx = torch.tensor(np.cos(np.outer(grid, m)))   # (n_grid, M+1)
        self.Cy = self.Cx.clone()
        c0 = self._fit_ic()
        self.c0 = c0

    def _fit_ic(self):
        X, Y = np.meshgrid(self.x_grid, self.y_grid, indexing="ij")
        vals = u0_2d(X, Y)
        m = np.arange(self.M + 1)
        Cx = np.cos(np.outer(self.x_grid, m))
        A = np.kron(Cx, Cx)
        c, *_ = np.linalg.lstsq(A, vals.flatten(), rcond=None)
        return torch.tensor(c.reshape(self.M + 1, self.M + 1))

    def field(self, C):
        """C: (B, M+1, M+1) -> u: (B, n_grid, n_grid)"""
        tmp = torch.einsum("xm,bmn->bxn", self.Cx, C)
        return torch.einsum("yn,bxn->bxy", self.Cy, tmp)

    def field_at(self, C, x_pts, y_pts):
        m = np.arange(self.M + 1)
        Cxp = torch.tensor(np.cos(np.outer(x_pts, m)))
        Cyp = torch.tensor(np.cos(np.outer(y_pts, m)))
        tmp = torch.einsum("xm,bmn->bxn", Cxp, C)
        return torch.einsum("yn,bxn->bxy", Cyp, tmp)


def train_fourier_spectral_pinn(eps=EPS2D, M=M_MODES, epochs=3000, nt_col=40,
                                 lr=2e-3, lam_ic=10.0, seed=42, warmup_epochs=500,
                                 n_grid=48, verbose=True):
    model = FourierSpectralPINN2D(eps, M=M, n_grid=n_grid, seed=seed)
    opt = torch.optim.Adam(model.net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    t_col = torch.linspace(1e-3, T2D, nt_col).reshape(-1, 1)
    t_ic = torch.zeros(1, 1)
    n_out = (M + 1) * (M + 1)

    warmup_opt = torch.optim.Adam(model.net.parameters(), lr=1e-3)
    t_warmup = torch.linspace(0, T2D, nt_col).reshape(-1, 1)
    c0_flat = model.c0.reshape(1, n_out)
    for _ in range(warmup_epochs):
        warmup_opt.zero_grad()
        C_pred = model.net(t_warmup).reshape(-1, n_out)
        loss_w = ((C_pred - c0_flat)**2).mean()
        loss_w.backward()
        warmup_opt.step()
    if verbose:
        print(f"  [Fourier-PINN warmup] final loss = {loss_w.item():.3e}")

    t0 = time.time()
    for ep in range(epochs):
        opt.zero_grad()
        t_col_ = t_col.clone().requires_grad_(True)
        C = model.net(t_col_)                       # (B,M+1,M+1)
        Cflat = C.reshape(-1, n_out)
        cols = []
        for k in range(n_out):
            gk = torch.autograd.grad(Cflat[:, k], t_col_, grad_outputs=torch.ones_like(Cflat[:, k]),
                                      create_graph=True, retain_graph=True)[0]
            cols.append(gk)
        dCdt = torch.cat(cols, dim=1).reshape(-1, M + 1, M + 1)

        u = model.field(C)
        u_t = model.field(dCdt)
        lap_C = model.lap_eig.unsqueeze(0) * C
        lap_u = model.field(lap_C)
        res = u_t - eps**2 * lap_u + f_ac(u)
        L_pde = (res**2).mean()

        C_ic = model.net(t_ic)
        L_ic = ((C_ic - model.c0)**2).mean()

        loss = L_pde + lam_ic * L_ic
        loss.backward()
        opt.step()
        sched.step()
        if verbose and (ep % 500 == 0 or ep == epochs - 1):
            print(f"  [Fourier-PINN] ep {ep:>5} loss={loss.item():.3e} pde={L_pde.item():.3e} ic={L_ic.item():.3e}")
    return model, time.time() - t0


def symmetry_violation_fourier(model, t_val, n=49):
    # NOTE: use an endpoint-INCLUSIVE grid so that grid[n-1-i] == -grid[i]
    # exactly; the model's *training* grid is periodic (endpoint=False),
    # but for this diagnostic we just need a grid where reflection is a
    # simple index reversal.
    grid = np.linspace(-np.pi, np.pi, n, endpoint=True)
    with torch.no_grad():
        C = model.net(torch.tensor([[t_val]]))
        u = model.field_at(C, grid, grid).numpy()[0]
    u_negx = u[::-1, :]
    u_negy = u[:, ::-1]
    S = np.mean((u - u_negx)**2) + np.mean((u - u_negy)**2)
    return S


# ------------------------------------------------------------------
# Vanilla 2D+t PINN with periodic (Fourier-feature) input encoding
# ------------------------------------------------------------------
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

    def forward(self, z):
        return self.net(z)


def encode(x, y, t):
    return torch.cat([torch.cos(x), torch.sin(x), torch.cos(y), torch.sin(y), t], dim=1)


def grad(y, x):
    return torch.autograd.grad(y, x, grad_outputs=torch.ones_like(y), create_graph=True)[0]


def train_vanilla_pinn_2d(eps=EPS2D, epochs=3000, n_xy=24, nt_col=15,
                           n_ic=30, lr=1e-3, lam_ic=10.0, seed=42, verbose=True):
    torch.manual_seed(seed)
    net = MLP([5, 64, 64, 64, 1])
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    xs = torch.linspace(-np.pi, np.pi, n_xy)
    ys = torch.linspace(-np.pi, np.pi, n_xy)
    ts = torch.linspace(1e-3, T2D, nt_col)
    Xc, Yc, Tc = torch.meshgrid(xs, ys, ts, indexing="ij")
    Xc = Xc.reshape(-1, 1).requires_grad_(True)
    Yc = Yc.reshape(-1, 1).requires_grad_(True)
    Tc = Tc.reshape(-1, 1).requires_grad_(True)

    xs_ic = torch.linspace(-np.pi, np.pi, n_ic)
    Xi, Yi = torch.meshgrid(xs_ic, xs_ic, indexing="ij")
    Xi = Xi.reshape(-1, 1)
    Yi = Yi.reshape(-1, 1)
    Ti = torch.zeros_like(Xi)
    Ui = torch.tensor(u0_2d(Xi.numpy(), Yi.numpy()))

    t0 = time.time()
    for ep in range(epochs):
        opt.zero_grad()
        z = encode(Xc, Yc, Tc)
        u = net(z)
        u_x = grad(u, Xc)
        u_xx = grad(u_x, Xc)
        u_y = grad(u, Yc)
        u_yy = grad(u_y, Yc)
        u_t = grad(u, Tc)
        res = u_t - eps**2 * (u_xx + u_yy) + f_ac(u)
        L_pde = (res**2).mean()

        z_ic = encode(Xi, Yi, Ti)
        u_ic_pred = net(z_ic)
        L_ic = ((u_ic_pred - Ui)**2).mean()

        loss = L_pde + lam_ic * L_ic
        loss.backward()
        opt.step()
        sched.step()
        if verbose and (ep % 500 == 0 or ep == epochs - 1):
            print(f"  [Vanilla-PINN-2D] ep {ep:>5} loss={loss.item():.3e} pde={L_pde.item():.3e} ic={L_ic.item():.3e}")
    return net, time.time() - t0


def symmetry_violation_vanilla(net, t_val, n=49):
    grid = np.linspace(-np.pi, np.pi, n, endpoint=True)
    X, Y = np.meshgrid(grid, grid, indexing="ij")
    with torch.no_grad():
        Xt = torch.tensor(X.flatten()).reshape(-1, 1)
        Yt = torch.tensor(Y.flatten()).reshape(-1, 1)
        Tt = torch.full_like(Xt, t_val)
        z = encode(Xt, Yt, Tt)
        u = net(z).numpy().reshape(n, n)
    u_negx = u[::-1, :]
    u_negy = u[:, ::-1]
    S = np.mean((u - u_negx)**2) + np.mean((u - u_negy)**2)
    return S


if __name__ == "__main__":
    print("=== Training Fourier-Spectral-PINN (2D, cos-cos basis) ===")
    fmodel, ftime = train_fourier_spectral_pinn(epochs=3000)
    print("=== Training vanilla PINN (2D, Fourier-feature input) ===")
    vnet, vtime = train_vanilla_pinn_2d(epochs=3000)

    print(f"\n{'t':>5} {'S_Fourier-Spectral-PINN':>26} {'S_Vanilla-PINN':>18}")
    for t_val in [0.0, 0.2, 0.5, 1.0]:
        Sf = symmetry_violation_fourier(fmodel, t_val)
        Sv = symmetry_violation_vanilla(vnet, t_val)
        print(f"{t_val:>5.1f} {Sf:>26.3e} {Sv:>18.3e}")
    print(f"\nTraining time: Fourier-Spectral-PINN={ftime:.0f}s, Vanilla-PINN-2D={vtime:.0f}s")
