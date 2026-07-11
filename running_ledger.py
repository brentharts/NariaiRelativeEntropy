"""
Running-ledger computation for Appendix C of
"Relative Entropy on the Nariai Horizon".

The cut-resolved (partial) modular flux
    R(u_c) = int_{u_c}^{+inf} t(u) du,
    t(u)   = kappa * (-U) <:T_UU:>_S  per du,
is the ledger accumulated between the cut u_c and the bifurcation
surface. Claims verified here:
  (C-a) R(-inf) = total = 2*sinh^2(r) * I1  (interference integrates to 0)
  (C-b) R(u_c) is non-monotonic: it dips inside negativity windows
  (C-c) every dip is bounded by the budget: max dip depth <= D(W)
Run:  python3 running_ledger.py
"""
import numpy as np

kappa = 1.0
w0, sg = 3.0, 0.7
w = np.linspace(1e-4, 12.0, 2400)
aw = np.exp(-(w - w0) ** 2 / (2 * sg ** 2))
u = np.linspace(-40.0, 40.0, 32001)
du = u[1] - u[0]
coef = (-1j * w * aw) * (w[1] - w[0])
dchi = np.empty(u.shape, dtype=complex)
for i0 in range(0, len(u), 4000):
    blk = u[i0:i0 + 4000]
    dchi[i0:i0 + 4000] = np.exp(-1j * np.outer(blk, w)) @ coef

I1 = np.trapezoid(np.abs(dchi) ** 2, dx=du) / kappa

print("cut-resolved ledger R(u_c);  total prediction = 2 sinh^2(r) I1")
for r, th in [(0.5, 0.0), (0.5, np.pi), (1.0, np.pi / 2)]:
    t = (2 * np.sinh(r) ** 2 * np.abs(dchi) ** 2
         - np.sinh(2 * r) * np.real(np.exp(-1j * th) * dchi ** 2)) / kappa
    # R(u_c): integrate from the right (bifurcation side) backward
    R = np.concatenate([[0.0], np.cumsum(t[::-1]) * du])[:-1][::-1]
    total = R[0]
    pred = 2 * np.sinh(r) ** 2 * I1
    # dips: places where R exceeds its running-from-the-right minimum envelope
    # depth of the largest dip = max over u_c of (max_{u>u_c} R  -  R(u_c))
    running_max_from_right = np.maximum.accumulate(R[::-1])[::-1]
    dip = np.max(running_max_from_right - R)
    D_W = -np.trapezoid(np.where(t < 0, t, 0.0), dx=du)
    print(f" r={r:.1f} th={th:5.2f}:  R(-inf) = {total: .6e}  "
          f"2sinh^2(r)I1 = {pred:.6e}  rel.diff = {abs(total-pred)/pred:.1e}")
    print(f"                largest dip = {dip:.6e}  D(W) = {D_W:.6e}  "
          f"dip <= D(W): {dip <= D_W * (1 + 1e-9)}")
