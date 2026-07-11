"""
Numerical verification for Appendix B of
"Relative Entropy on the Nariai Horizon".

Three checks:
  (B1) For boost-positive-frequency horizon packets, the interference
       functional I2 vanishes (exactly, in the continuum), so the
       squeeze-shape inequality |I2|/I1 <= tanh(r) holds trivially.
       Mechanism: with U = -(1/k) e^{-k u}, the modular weight obeys
       (-U) dU = (1/k) du, so I2 = (1/k) * [Fourier transform of
       (d_u chi)^2 at total frequency zero] = 0 for packets built
       from strictly positive boost frequencies.
  (B2) The squeezed-state flux <T_UU> has genuinely negative local
       windows, while the modular-weighted total (= S_rel / 2pi)
       stays positive: the anti-evaporation budget in action.
  (B3) One-mode Gaussian check of the sum rule S_rel = beta*Delta E
       (i.e. Delta S = 0 for wedge-local unitary excitations),
       via covariance-matrix formulas.

kappa = 1 throughout (Nariai: kappa = sqrt(Lambda)).
"""
import numpy as np

kappa = 1.0

# ----------------------------------------------------------------------
# (B1) Boost-positive-frequency packet: I1, I2 in boost coordinates
# ----------------------------------------------------------------------
# chi(u) = int_0^inf dw a(w) e^{-i w u},  a(w) Gaussian at w0, width sg
w0, sg = 3.0, 0.7
w = np.linspace(1e-4, 12.0, 2400)          # strictly positive frequencies
aw = np.exp(-(w - w0) ** 2 / (2 * sg ** 2))
u = np.linspace(-40.0, 40.0, 32001)
du = u[1] - u[0]

coef = (-1j * w * aw) * (w[1] - w[0])
dchi_du = np.empty(u.shape, dtype=complex)
for i0 in range(0, len(u), 4000):                 # chunked to bound memory
    blk = u[i0:i0 + 4000]
    dchi_du[i0:i0 + 4000] = np.exp(-1j * np.outer(blk, w)) @ coef

I1 = (1.0 / kappa) * np.trapezoid(np.abs(dchi_du) ** 2, dx=du)
I2 = (1.0 / kappa) * np.trapezoid(dchi_du ** 2, dx=du)
ratio_pos = abs(I2) / I1

# Cross-check in affine U coordinates: U = -(1/k) e^{-k u}
# (-U)|d_U chi|^2 dU = (1/k)|d_u chi|^2 du  -- verify numerically
U = -(1.0 / kappa) * np.exp(-kappa * u)
dchi_dU = dchi_du * np.exp(kappa * u)      # d_U = e^{k u} d_u  (dU/du = e^{-k u})
I1_U = np.trapezoid((-U) * np.abs(dchi_dU) ** 2 * np.gradient(U), dx=1.0)  # dU = grad
# cleaner: integrate in u with Jacobian
I1_U = np.trapezoid((-U) * np.abs(dchi_dU) ** 2 * np.exp(-kappa * u), dx=du)
I2_U = np.trapezoid((-U) * (dchi_dU ** 2) * np.exp(-kappa * u), dx=du)

print("== (B1) boost-positive-frequency packet ==")
print(f"I1 (boost coords)          = {I1:.6e}")
print(f"I1 (affine U cross-check)  = {I1_U:.6e}   rel.diff = {abs(I1-I1_U)/I1:.2e}")
print(f"|I2|/I1                    = {ratio_pos:.3e}   (continuum prediction: 0)")
print(f"|I2_U|/I1                  = {abs(I2_U)/I1:.3e}")

# Contrast: a packet NOT adapted to the wedge (real Gaussian in U).
# Real data => (d_U chi)^2 = |d_U chi|^2 pointwise => |I2|/I1 = 1 exactly.
U0, s = -2.0, 0.5
Ugrid = np.linspace(-8.0, -1e-3, 200001)
chiR = np.exp(-(Ugrid - U0) ** 2 / (2 * s ** 2))
dchiR = np.gradient(chiR, Ugrid)
I1R = np.trapezoid((-Ugrid) * dchiR ** 2, Ugrid)
I2R = np.trapezoid((-Ugrid) * dchiR ** 2, Ugrid)   # real => identical
print("\n== (B1') real (non-boost-adapted) packet ==")
print(f"|I2|/I1 = {abs(I2R)/I1R:.6f}   => tanh(r) >= 1 required: no wedge-local squeeze exists")

# ----------------------------------------------------------------------
# (B2) Negativity windows vs. positive modular-weighted total
# ----------------------------------------------------------------------
print("\n== (B2) squeezed flux: local negativity, global positivity ==")
for r, th in [(0.5, 0.0), (0.5, np.pi), (1.0, np.pi / 2)]:
    # modular-weighted flux density in boost time:  t(u) du = 2pi-normalized
    t = (2 * np.sinh(r) ** 2 * np.abs(dchi_du) ** 2
         - np.sinh(2 * r) * np.real(np.exp(-1j * th) * dchi_du ** 2)) / kappa
    total = np.trapezoid(t, dx=du)                    # = kappa-normalized S_rel/2pi
    neg = -np.trapezoid(np.where(t < 0, t, 0.0), dx=du)   # funded depth D(W)
    pos = np.trapezoid(np.where(t > 0, t, 0.0), dx=du)
    frac_neg_support = np.mean(t < 0)
    print(f" r={r:.1f} th={th:5.2f}:  S_rel/2pi = {total: .5e}  "
          f"D(W) = {neg:.5e}  D(W)/pos-part = {neg/pos:.3f}  "
          f"neg support = {100*frac_neg_support:.1f}%  total>0: {total > 0}")

# ----------------------------------------------------------------------
# (B3) One-mode Gaussian check of the sum rule  S_rel = beta * Delta E
# ----------------------------------------------------------------------
# Wedge restriction of the global vacuum = thermal state at beta = 2pi/kappa
# for a boost mode of frequency w_m. Wedge-local squeeze = unitary on that
# mode. Covariance convention: sigma_vac = I, <n> = (Tr sigma - 2)/4.
print("\n== (B3) one-mode sum-rule check ==")
def SvN(nu):
    a, b = (nu + 1) / 2, (nu - 1) / 2
    return a * np.log(a) - (b * np.log(b) if b > 0 else 0.0)

beta = 2 * np.pi / kappa
for w_m, r in [(0.5, 0.3), (1.0, 0.8), (2.0, 1.5)]:
    nu = 1.0 / np.tanh(beta * w_m / 2)             # thermal symplectic eigenvalue
    sig_th = nu * np.eye(2)
    S = np.diag([np.exp(r), np.exp(-r)])           # squeeze, theta = 0
    sig_sq = S @ sig_th @ S.T
    # symplectic eigenvalue of sig_sq (should equal nu: unitary invariance)
    nu_sq = np.sqrt(np.linalg.det(sig_sq))
    E_th = w_m * np.trace(sig_th) / 4
    E_sq = w_m * np.trace(sig_sq) / 4
    dE = E_sq - E_th
    dS_vN = SvN(nu_sq) - SvN(nu)                   # should be 0
    # Araki/Umegaki relative entropy vs thermal reference:
    lnZ = -np.log(2 * np.sinh(beta * w_m / 2))
    S_rel = beta * E_sq + lnZ - SvN(nu_sq)
    S_th_check = beta * E_th + lnZ - SvN(nu)       # = 0 identically
    print(f" w={w_m:.1f} r={r:.1f}: nu_sq-nu = {nu_sq-nu:.2e}  dS_vN = {dS_vN:.2e}  "
          f"S_rel = {S_rel:.6f}  beta*dE = {beta*dE:.6f}  "
          f"|S_rel - beta*dE| = {abs(S_rel-beta*dE):.2e}")
    # analytic form: beta * w * sinh^2(r) * coth(beta w / 2) ... times 2? check:
    analytic = beta * w_m * nu * np.sinh(r) ** 2
    print(f"           analytic beta*w*nu*sinh^2(r) = {analytic:.6f}")
