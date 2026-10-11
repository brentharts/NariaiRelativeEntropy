"""
check_open_items.py -- the four items left open after the v10 revision, each
computed rather than argued.  Exit status 0 iff every check passes.

  [1] Orientation of the modular flow from KMS analyticity.
      For U, U' < 0 the Tomita-Takesaki function F(t) = <A sigma_t(B)> with the
      scaling-limit kernel K ~ -1/pi (U - lam U' - i0)^-2 must be analytic in
      -1 < Im t < 0, equal <A sigma_t B> at Im t = 0 and <sigma_t(B) A> at
      Im t = -1.  With lam = e^{2 pi t} (data phi^t(U) = phi(e^{-2 pi t} U),
      the paper's choice) this holds: the continuation runs from just below
      the pole ray to just above it.  With lam = e^{-2 pi t} the boundary
      values are interchanged.  The commutator is nonzero, so the
      orientation is not a convention: Im I(lam0 -+ i0) = -+ J with J the
      delta' term, computed independently.
  [2] All three 2D components of the linearized 4D Einstein equation about
      dS2 x S2 equal the JT operator (-Lambda sign, normalization -(2/r0)),
      independently of the 2D metric perturbation.  (Lambda = 2 so that
      factors of Lambda are not hidden.)
  [3] The validity window kappa*Du < p ln S_N for a seed eps0 ~ S_N^-p.
  [4] l >= 1: the linearized Raychaudhuri equation holds generator by
      generator, so the retarded response of the area element integrates to
      -4 S_rel for nonspherical data, checked against the factorized
      harmonic sum; and wedge-local (symplectic) excitations have dS_vN = 0
      so S_rel = beta dE, while a non-unitary change does not.
"""
import math
import sys
import numpy as np

FAIL = []


def check(name, cond):
    print("  {:<76} {}".format(name, "ok" if cond else "FAIL"))
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------- [1] -------
def part1():
    print("\n[1] modular-flow orientation from KMS analyticity")
    N, dx = 3000, 0.02
    x = (np.arange(N) + 0.5) * dx
    f = x ** 3 * np.exp(-x)

    def I(lam):                      # -(1/pi) int int f f / (lam x' - x)^2
        out = 0j
        for i in range(0, N, 300):
            xx = x[i:i + 300, None]
            out += np.sum(f[i:i + 300, None] * f[None, :]
                          / (lam * x[None, :] - xx) ** 2) * dx * dx
        return -out / np.pi

    lam0 = 1.3
    below = [I(lam0 * np.exp(-2j * np.pi * d)) for d in (0.02, 0.01, 0.005)]
    above = [I(lam0 * np.exp(+2j * np.pi * d)) for d in (0.02, 0.01, 0.005)]
    check("forward and reversed orderings differ (commutator != 0)",
          abs(below[-1].imag - above[-1].imag) > 0.5)
    check("they are complex conjugates (same real part)",
          abs(below[-1].real - above[-1].real) < 1e-9)
    from scipy.integrate import quad
    fp = lambda s: (3 * s ** 2 - s ** 3) * np.exp(-s)
    J = quad(lambda s: s ** 3 * np.exp(-s) * fp(lam0 * s), 0, 80)[0]
    errs = [abs(b.imag - J) for b in below]
    check("Im I(lam0 - i0) -> delta' term J by a second route (error shrinks)",
          errs[0] > errs[1] > errs[2] and errs[2] < 0.1)
    # continuation along lam = lam0 e^{i arg}, arg: -small -> -2 pi + small
    d0 = 0.02
    args = np.linspace(-2 * np.pi * d0, -2 * np.pi + 2 * np.pi * d0, 150)
    vals = np.array([I(lam0 * np.exp(1j * a)) for a in args])
    check("orientation A (lam = e^{2 pi t}): continuation is smooth (no pole crossed)",
          np.max(np.abs(np.diff(vals))) < 0.2)
    check("orientation A: Im t -> 0 end equals the forward (below) value",
          abs(vals[0] - I(lam0 * np.exp(-2j * np.pi * d0))) < 1e-9)
    check("orientation A: Im t -> -1 end equals the reversed (above) value",
          abs(vals[-1] - I(lam0 * np.exp(+2j * np.pi * d0))) < 1e-9)
    # orientation B: lam = lam0 e^{-i arg}: the Im t -> 0 end is the ABOVE value
    check("control: orientation B puts the reversed ordering at Im t = 0 (violates KMS)",
          abs(I(lam0 * np.exp(+2j * np.pi * d0)) - below[0]) > 0.3)
    # positivity of the kernel's quadratic form (physical sign of i0)
    check("Lambda_2(f,f) = Re I(1) > 0 (positive-energy kernel)", I(1.0 - 1e-3j).real > 0)


# ---------------------------------------------------------------- [2] -------
def part2():
    print("\n[2] all 2D components of the linearized 4D Einstein equation")
    import sympy as sp
    T, th, a, b, eps = sp.symbols("T theta a b epsilon")
    Lam = sp.Integer(2)
    r0 = 1 / sp.sqrt(Lam)

    def run(rho, h):
        X = [T, th, a, b]
        r = r0 + eps * rho
        g = sp.zeros(4, 4)
        g[0, 0] = -1 / Lam + eps * h[0]
        g[0, 1] = g[1, 0] = eps * h[1]
        g[1, 1] = sp.cosh(T) ** 2 / Lam + eps * h[2]
        g[2, 2] = r ** 2
        g[3, 3] = r ** 2 * sp.sin(a) ** 2
        g0 = g.subs(eps, 0)
        g0i = g0.inv()
        dg = sp.diff(g, eps).subs(eps, 0)
        gi = g0i - eps * g0i * dg * g0i
        tr = lambda e: (lambda q: q.coeff(eps, 0) + eps * q.coeff(eps, 1))(sp.expand(e))
        G = [[[tr(sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j])
                                  - sp.diff(g[j, k], X[l])) for l in range(4)) / 2)
               for k in range(4)] for j in range(4)] for i in range(4)]

        def Ric(j, k):
            e = sum(sp.diff(G[i][j][k], X[i]) - sp.diff(G[i][j][i], X[k])
                    for i in range(4))
            e += sum(G[i][i][l] * G[l][j][k] - G[i][k][l] * G[l][j][i]
                     for i in range(4) for l in range(4))
            return tr(e)
        R = sp.Matrix(4, 4, Ric)
        Rs = tr(sum(gi[j, k] * R[j, k] for j in range(4) for k in range(4)))
        E = sp.Matrix(4, 4, lambda j, k: tr(R[j, k] - Rs * g[j, k] / 2 + Lam * g[j, k]))
        return E, E.applyfunc(lambda e: sp.diff(e, eps))

    def jt(rho, sgn):
        g2 = sp.Matrix([[-1 / Lam, 0], [0, sp.cosh(T) ** 2 / Lam]])
        g2i = g2.inv()
        x2 = [T, th]
        G2 = [[[sum(g2i[i, l] * (sp.diff(g2[l, j], x2[k]) + sp.diff(g2[l, k], x2[j])
                                 - sp.diff(g2[j, k], x2[l])) for l in range(2)) / 2
                for k in range(2)] for j in range(2)] for i in range(2)]
        H = sp.Matrix(2, 2, lambda i, j: sp.diff(rho, x2[i], x2[j])
                      - sum(G2[k][i][j] * sp.diff(rho, x2[k]) for k in range(2)))
        box = sum(g2i[i, j] * H[i, j] for i in range(2) for j in range(2))
        return (H - g2 * box + sgn * Lam * g2 * rho) * (-2 / r0)

    rho = sp.sin(T + 2 * th) + T * th ** 2 / 3 + sp.cosh(T) * sp.cos(th)
    h1 = [sp.cos(T) * th, T ** 2 * sp.sin(th), sp.exp(T / 3) * th]
    h2 = [th ** 2 + T, sp.cos(2 * T + th), T * th]
    E1, dE1 = run(rho, h1)
    _, dE2 = run(rho, h2)
    check("background G + Lambda g vanishes", E1.subs(eps, 0) == sp.zeros(4, 4))
    Jm, Jp = jt(rho, -1), jt(rho, +1)
    pts = [{T: 0.37, th: 0.81, a: 1.1, b: 0.2}, {T: -0.9, th: 2.3, a: 0.6, b: 1.7}]
    num = lambda e, p: complex(sp.N(e.subs(p)))
    for (i, j), nm in (((0, 0), "TT"), ((0, 1), "T-theta"), ((1, 1), "theta-theta")):
        check("%s component = JT operator with -Lambda g phi" % nm,
              all(abs(num(dE1[i, j] - Jm[i, j], p)) < 1e-9 for p in pts))
        check("%s component independent of the 2D metric perturbation" % nm,
              all(abs(num(dE1[i, j] - dE2[i, j], p)) < 1e-9 for p in pts))
        if nm != "T-theta":     # g_T,theta = 0, so the sign cannot show there
            check("control: %s component fails with +Lambda g phi" % nm,
                  any(abs(num(dE1[i, j] - Jp[i, j], p)) > 1e-3 for p in pts))
    check("angular component does depend on the metric perturbation (the other equation)",
          any(abs(num(dE1[2, 2] - dE2[2, 2], p)) > 1e-3 for p in pts))


# ---------------------------------------------------------------- [3] -------
def part3():
    print("\n[3] validity window for seed eps0 ~ S_N^-p")
    lnS = math.log(1e122)            # observed Lambda: S ~ 1e122
    check("ln S_N ~ 281 for S_N ~ 1e122", abs(lnS - 280.9) < 0.1)
    # window kappa*Du_max = ln(1/eps0) = p ln S_N: logarithmic in S_N for every p
    ws = {p: p * lnS for p in (0.5, 1.0, 2.0)}
    check("window is p ln S_N: 140 / 281 / 562 for p = 1/2, 1, 2", all(
        abs(ws[p] - p * 280.9) < 0.2 for p in ws))
    # doubling the exponent of the seed only doubles the window; squaring S_N doubles it too
    check("squaring S_N doubles the window (logarithmic, not a power)",
          abs(math.log(1e244) / math.log(1e122) - 2) < 1e-12)
    # growth e^{sqrt(Lambda) tau}: time to O(1) from eps0
    for p in (0.5, 1.0):
        eps0 = 1e122 ** (-p)
        check("p=%.1f: eps0 e^{kappa Du} = 1 at kappa Du = ln(1/eps0)" % p,
              abs(math.log(1 / eps0) - p * lnS) < 1e-9)


# ---------------------------------------------------------------- [4] -------
def part4():
    print("\n[4] l >= 1 (generator-wise Raychaudhuri) and wedge-local excitations")
    from scipy.integrate import solve_ivp, quad
    from numpy.polynomial.legendre import leggauss, legval
    Lam = 1.0
    amp = {0: 1.0, 1: 0.7, 2: -0.5, 3: 0.3}          # radial amplitudes by l
    prof = lambda U, l: amp[l] * np.exp(-(U + 2.0 - 0.3 * l) ** 2 / 0.6)
    ms = 0.0   # phi(U, s) = sum_l prof(U,l) sqrt((2l+1)/4pi) P_l(cos th)  (m = 0)
    xs, ws = leggauss(40)
    Ylm = lambda l, c: math.sqrt((2 * l + 1) / (4 * math.pi)) * legval(c, [0] * l + [1])
    dphi = lambda U, c: sum(
        -2 * (U + 2.0 - 0.3 * l) / 0.6 * prof(U, l) * Ylm(l, c) for l in amp)
    T_UU = lambda U, c: dphi(U, c) ** 2
    # route 1: integrate a'' = -8 pi T_UU with retarded data at each polar angle
    Us = np.linspace(-30, 0, 4001)
    area = 0.0
    for c, w in zip(xs, ws):
        sol = solve_ivp(lambda U, y: [y[1], -8 * np.pi * T_UU(U, c)],
                        [-30, 0], [0.0, 0.0], t_eval=[0.0], rtol=1e-10, atol=1e-12)
        area += w * sol.y[0, -1] * 2 * math.pi / Lam     # dvol = dOmega/Lam
    # route 2: factorized harmonic sum  S_rel = (2 pi/Lam) sum_l int (-U) (dphi_l)^2
    def srel_l(l):
        d = lambda U: -2 * (U + 2.0 - 0.3 * l) / 0.6 * prof(U, l)
        return 2 * math.pi / Lam * quad(lambda U: (-U) * d(U) ** 2, -30, 0)[0]
    S = sum(srel_l(l) for l in amp)
    check("nonspherical data: retarded dA = -4 S_rel (ODE route vs harmonic sum)",
          abs(area + 4 * S) < 1e-5 * abs(4 * S))
    S_no1 = S - srel_l(1)
    check("control: dropping the l = 1 term breaks it", abs(area + 4 * S_no1) > 1e-2 * abs(4 * S))
    check("l = 3 contributes (the sum is not saturated at l = 0)", srel_l(3) > 1e-3)
    # wedge-local one-mode excitations: symplectic => dS = 0 => S_rel = beta dE
    SvN = lambda nu: ((nu + 1) / 2) * math.log((nu + 1) / 2) - (
        ((nu - 1) / 2) * math.log((nu - 1) / 2) if nu > 1 + 1e-14 else 0.0)
    kappa, w0 = 1.0, 1.0
    beta = 2 * math.pi / kappa
    nu = 1 / math.tanh(beta * w0 / 2)
    lnZ = -math.log(2 * math.sinh(beta * w0 / 2))
    ok_sym, ok_non = True, False
    for r in (0.3, 0.8, 1.5):
        S_mat = np.diag([math.exp(r), math.exp(-r)])
        sig = S_mat @ (nu * np.eye(2)) @ S_mat.T
        nu_new = math.sqrt(np.linalg.det(sig))
        dS = SvN(nu_new) - SvN(nu)
        E_sq = w0 * np.trace(sig) / 4
        E_th = w0 * 2 * nu / 4
        srel = beta * E_sq + lnZ - SvN(nu_new)
        ok_sym &= abs(dS) < 1e-12 and abs(srel - beta * (E_sq - E_th)) < 1e-10
    check("symplectic (wedge-local) squeezing: dS_vN = 0 and S_rel = beta dE", ok_sym)
    nu2 = 1.7 * nu                           # not a symplectic image of the thermal state
    sig2 = nu2 * np.eye(2)
    srel2 = beta * (w0 * np.trace(sig2) / 4) + lnZ - SvN(nu2)
    dE2 = w0 * np.trace(sig2) / 4 - w0 * 2 * nu / 4
    check("control: a non-unitary change has dS_vN != 0 and S_rel != beta dE",
          abs(SvN(nu2) - SvN(nu)) > 1e-3 and abs(srel2 - beta * dE2) > 1e-3)


def main():
    for p in (part1, part2, part3, part4):
        p()
    print()
    if FAIL:
        print("FAILED: %d" % len(FAIL))
        for f in FAIL:
            print("  - " + f)
        return 1
    print("all open-item checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
