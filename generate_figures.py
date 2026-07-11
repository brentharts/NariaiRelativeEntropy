"""
Generate figures for the supplementary material of
"Relative Entropy on the Nariai Horizon".

Outputs (into figures/):
  fig_flux_windows.pdf : squeezed horizon flux in boost time, with the
                         negativity windows shaded -- local NEC violation
                         funded by positive flux elsewhere (zero-sum).
  fig_dichotomy.pdf    : the squeeze-shape functional |I2|/I1 as the packet's
                         spectral weight leaks to negative boost frequencies;
                         the exact wedge-locality marker.
  fig_one_mode.pdf     : one-mode sum-rule check, S_rel = beta*Delta E,
                         covariance-matrix points vs the closed form
                         beta*w*sinh^2(r)*coth(beta*w/2).
Run:  python3 generate_figures.py
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("figures", exist_ok=True)
kappa = 1.0

# ---------------------------------------------------------------- packet ---
def packet_dchi(u, w0=3.0, sg=0.7, wmax=12.0, nw=2400):
    """d_u chi for a boost-positive-frequency Gaussian spectral packet."""
    w = np.linspace(1e-4, wmax, nw)
    aw = np.exp(-(w - w0) ** 2 / (2 * sg ** 2))
    coef = (-1j * w * aw) * (w[1] - w[0])
    out = np.empty(u.shape, dtype=complex)
    for i0 in range(0, len(u), 4000):
        blk = u[i0:i0 + 4000]
        out[i0:i0 + 4000] = np.exp(-1j * np.outer(blk, w)) @ coef
    return out

# ------------------------------------------------- fig 1: flux windows ----
u = np.linspace(-4.0, 4.0, 8001)
dchi = packet_dchi(u)
r, th = 0.5, np.pi
t = (2 * np.sinh(r) ** 2 * np.abs(dchi) ** 2
     - np.sinh(2 * r) * np.real(np.exp(-1j * th) * dchi ** 2)) / kappa

fig, ax = plt.subplots(figsize=(7.0, 3.2))
ax.plot(u, t, lw=1.2, color="k")
ax.fill_between(u, t, 0, where=(t < 0), color="crimson", alpha=0.45,
                label=r"negativity windows $\mathcal{W}$")
ax.fill_between(u, t, 0, where=(t > 0), color="steelblue", alpha=0.25,
                label="positive (prepaying) flux")
ax.axhline(0, color="gray", lw=0.6)
ax.set_xlabel(r"boost time $u$")
ax.set_ylabel(r"$\kappa\,(-U)\,\langle{:}T_{UU}{:}\rangle_S\;$ (per $du$)")
ax.set_title(rf"Squeezed horizon flux, $r={r}$, $\theta=\pi$: "
             r"local negativity, zero-sum interference")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("figures/fig_flux_windows.pdf")
plt.close(fig)

# --------------------------------------------- fig 2: locality marker -----
# Let spectral weight leak below w=0: a(w) Gaussian centered w0, we sweep
# w0/sg downward; the negative-frequency content grows and |I2|/I1 rises
# from ~0 toward 1 (real-packet limit).
uu = np.linspace(-60.0, 60.0, 32001)
duu = uu[1] - uu[0]
ratios, leaks = [], []
sg = 1.0
for w0 in np.linspace(3.0, -1.5, 25):
    w = np.linspace(-8.0, 12.0, 4000)            # allow negative frequencies
    aw = np.exp(-(w - w0) ** 2 / (2 * sg ** 2))
    coef = (-1j * w * aw) * (w[1] - w[0])
    d = np.empty(uu.shape, dtype=complex)
    for i0 in range(0, len(uu), 4000):
        blk = uu[i0:i0 + 4000]
        d[i0:i0 + 4000] = np.exp(-1j * np.outer(blk, w)) @ coef
    I1 = np.trapezoid(np.abs(d) ** 2, dx=duu)
    I2 = np.trapezoid(d ** 2, dx=duu)
    ratios.append(abs(I2) / I1)
    leaks.append(np.trapezoid(aw[w < 0] ** 2, dx=w[1]-w[0])
                 / np.trapezoid(aw ** 2, dx=w[1]-w[0]))

fig, ax = plt.subplots(figsize=(6.4, 3.4))
ax.plot(leaks, ratios, "o-", ms=4, color="k")
ax.set_xlabel(r"negative-boost-frequency fraction of $|a(\omega)|^2$")
ax.set_ylabel(r"$|I_2|/I_1$")
ax.set_title("The squeeze--shape functional as an exact wedge-locality marker")
ax.axhline(0, color="steelblue", lw=0.8, ls="--")
ax.axhline(1, color="crimson", lw=0.8, ls="--")
ax.text(0.02, 0.06, "boost-adapted: $I_2=0$", fontsize=9, color="steelblue")
ax.text(0.30, 0.90, "boost-blind limit: $|I_2|/I_1=1$", fontsize=9, color="crimson")
fig.tight_layout()
fig.savefig("figures/fig_dichotomy.pdf")
plt.close(fig)

# ------------------------------------------------ fig 3: one-mode law -----
def SvN(nu):
    a, b = (nu + 1) / 2, (nu - 1) / 2
    return a * np.log(a) - (b * np.log(b) if b > 1e-14 else 0.0)

beta = 2 * np.pi / kappa
rr = np.linspace(0, 1.6, 33)
fig, ax = plt.subplots(figsize=(6.4, 3.4))
for w_m, c in [(0.5, "steelblue"), (1.0, "k"), (2.0, "crimson")]:
    nu = 1.0 / np.tanh(beta * w_m / 2)
    pts = []
    for r_ in rr:
        S = np.diag([np.exp(r_), np.exp(-r_)])
        sig = S @ (nu * np.eye(2)) @ S.T
        E_sq = w_m * np.trace(sig) / 4
        lnZ = -np.log(2 * np.sinh(beta * w_m / 2))
        pts.append(beta * E_sq + lnZ - SvN(np.sqrt(np.linalg.det(sig))))
    ax.plot(rr, beta * w_m * nu * np.sinh(rr) ** 2, "-", color=c, lw=1.0,
            label=rf"$\omega={w_m}$ (closed form)")
    ax.plot(rr[::4], np.array(pts)[::4], "o", color=c, ms=4)
ax.set_xlabel(r"squeezing parameter $r$")
ax.set_ylabel(r"$S^{\mathrm{rel}}$")
ax.set_title(r"One-mode sum rule: covariance points vs "
             r"$\beta\,\omega\,\sinh^2\!r\,\coth(\beta\omega/2)$")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("figures/fig_one_mode.pdf")
plt.close(fig)

print("wrote figures/fig_flux_windows.pdf, fig_dichotomy.pdf, fig_one_mode.pdf")
