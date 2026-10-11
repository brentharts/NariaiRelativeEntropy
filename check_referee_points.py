"""
check_referee_points.py -- the checkable claims of the referee report on v9,
each computed rather than argued.  Exit status 0 iff every check passes.

  [1] dS2 dilaton equation:  with g = Lambda^-1(-dT^2 + cosh^2 T dtheta^2)
      (R = 2 Lambda) the equation  nabla_a nabla_b phi - g_ab box phi
      + s Lambda g_ab phi = 0  is solved by sinh T, cosh T cos(theta),
      cosh T sin(theta) for s = -1 and for NO such mode when s = +1.
      These modes grow like e^{sqrt(Lambda) tau} in proper time tau: the
      growth rate is kappa, not kappa/S_N.
  [2] Sign of the area response.  phi'' = -(8 pi/Lambda) T_UU, T_UU >= 0.
      Retarded data (phi, phi' -> 0 at U -> -inf): phi(0) <= 0, i.e.
      dA = -4 S_rel (focusing).  Teleological data (phi'(0) = 0): the
      modular-weighted term still enters with a NEGATIVE sign, accompanied
      by a positive memory term proportional to (-U0) * (absorbed flux).
  [3] SdS horizon areas at fixed Lambda: the cosmological horizon area
      decreases monotonically from 12 pi/Lambda (M = 0) to 4 pi/Lambda
      (Nariai); the black-hole horizon area increases from 0 to 4 pi/Lambda.
      Nariai is the maximal black-hole and minimal cosmological horizon, and
      pure de Sitter has the largest entropy 3 pi/Lambda in the family.
"""
import sys
import numpy as np
import sympy as sp
from scipy.integrate import quad

FAIL = []


def check(name, cond):
    print("  {:<72} {}".format(name, "ok" if cond else "FAIL"))
    if not cond:
        FAIL.append(name)


def part1():
    print("\n[1] dS2 dilaton equation, sign of the Lambda term")
    T, th, L = sp.symbols("T theta Lambda", positive=True)
    x = [T, th]
    g = sp.Matrix([[-1 / L, 0], [0, sp.cosh(T) ** 2 / L]])
    gi = g.inv()
    G = [[[sum(gi[a, d] * (sp.diff(g[d, b], x[c]) + sp.diff(g[d, c], x[b])
                           - sp.diff(g[b, c], x[d])) for d in range(2)) / 2
           for c in range(2)] for b in range(2)] for a in range(2)]

    def hess(f):
        return sp.Matrix(2, 2, lambda a, b: sp.diff(f, x[a], x[b])
                         - sum(G[c][a][b] * sp.diff(f, x[c]) for c in range(2)))

    def box(f):
        H = hess(f)
        return sp.simplify(sum(gi[a, b] * H[a, b]
                               for a in range(2) for b in range(2)))

    def riem(a, b, c, d):
        return (sp.diff(G[a][b][d], x[c]) - sp.diff(G[a][b][c], x[d])
                + sum(G[a][c][e] * G[e][b][d] - G[a][d][e] * G[e][b][c]
                      for e in range(2)))
    ric = sp.Matrix(2, 2, lambda b, d: sum(riem(a, b, a, d) for a in range(2)))
    R = sp.simplify(sum(gi[b, d] * ric[b, d]
                        for b in range(2) for d in range(2)))
    check("R = 2 Lambda", sp.simplify(R - 2 * L) == 0)

    def solves(f, s):
        E = (hess(f) - g * box(f) + s * L * g * f).applyfunc(sp.simplify)
        return E == sp.zeros(2, 2)
    modes = [sp.sinh(T), sp.cosh(T) * sp.cos(th), sp.cosh(T) * sp.sin(th)]
    for i, f in enumerate(modes):
        check("mode %d solves the equation with -Lambda g phi" % i, solves(f, -1))
        check("control: mode %d fails with +Lambda g phi" % i, not solves(f, +1))
    check("control: tanh T is not a solution (either sign)",
          not solves(sp.tanh(T), -1) and not solves(sp.tanh(T), 1))
    # growth: phi = sinh T, proper time tau = T/sqrt(Lambda)
    check("sinh T grows like exp(sqrt(Lambda) tau): d ln/dtau -> sqrt(Lambda)",
          sp.limit(sp.diff(sp.log(sp.sinh(T)), T), T, sp.oo) == 1)


def part2():
    print("\n[2] sign of the area response")
    Lam = 1.0
    flux = lambda U: np.exp(-(U + 2.0) ** 2 / 0.5)
    src = lambda U: -(8 * np.pi / Lam) * flux(U)
    mod = quad(lambda U: (-U) * flux(U), -40, 0)[0]
    ret = quad(lambda U: (-U) * src(U), -40, 0)[0]
    check("retarded: phi(0) = -(8 pi/Lambda) * modular integral",
          abs(ret + 8 * np.pi / Lam * mod) < 1e-9)
    check("retarded: phi(0) < 0 (focusing, dA = -4 S_rel)", ret < 0)
    for U0 in (-3.0, -6.0, -12.0):
        d = -quad(lambda Up: src(Up) * (Up - U0), U0, 0)[0]
        mem = (8 * np.pi / Lam) * (-U0) * quad(flux, U0, 0)[0]
        mw = -(8 * np.pi / Lam) * quad(lambda Up: (-Up) * flux(Up), U0, 0)[0]
        check("teleological U0=%g: phi(0)-phi(U0) = memory + modular" % U0,
              abs(d - (mem + mw)) < 1e-8)
        check("teleological U0=%g: modular term negative, memory positive" % U0,
              mw < 0 < mem)
    check("control: the + sign of the paper's Eq. (48) is not the retarded one",
          np.sign(ret) != np.sign(+8 * np.pi / Lam * mod))


def part3():
    print("\n[3] Schwarzschild-de Sitter horizon areas at fixed Lambda = 1")
    Ms = np.linspace(1e-4, 1 / 3 - 1e-9, 1500)
    rb, rc = [], []
    for M in Ms:
        r = np.roots([-1 / 3, 0, 1, -2 * M])
        r = np.sort(r[np.isreal(r)].real)
        r = r[r > 0]
        rb.append(r[0])
        rc.append(r[-1])
    rb, rc = np.array(rb), np.array(rc)
    check("cosmological area -> 12 pi/Lambda as M -> 0",
          abs(4 * np.pi * rc[0] ** 2 - 12 * np.pi) < 0.1)
    check("cosmological area -> 4 pi/Lambda at Nariai",
          abs(4 * np.pi * rc[-1] ** 2 - 4 * np.pi) < 5e-3)   # double root: error ~ sqrt(1e-9)
    check("cosmological area strictly decreasing in M", bool(np.all(np.diff(rc) < 0)))
    check("black-hole area strictly increasing, 0 -> 4 pi/Lambda",
          bool(np.all(np.diff(rb) > 0)) and abs(4 * np.pi * rb[-1] ** 2 - 4 * np.pi) < 5e-3)
    check("pure de Sitter entropy 3 pi/Lambda exceeds Nariai pi/Lambda", 3 * np.pi > np.pi)
    check("control: 'no horizon exceeds the Nariai area' is false",
          4 * np.pi * rc[0] ** 2 > 4 * np.pi * rc[-1] ** 2)


def main():
    for p in (part1, part2, part3):
        p()
    print()
    if FAIL:
        print("FAILED: %d" % len(FAIL))
        return 1
    print("all referee-point checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
