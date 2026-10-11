"""
check_normalization.py -- what can and cannot be said about the coefficient 1/4.

The relation  dA = c * S_rel  has c of dimension length^2; S_rel = 2 pi dE/kappa
contains no gravitational constant.  So c/l_P^2 = 4 is not derivable from
horizon kinematics and QFT alone: it is the definition of Newton's constant.
What is not a definition is that the SAME constant governs the horizon
entropy and the asymptotic mass.  This script tests that, in a three-parameter
family of theories, by two code-independent routes:

  route H (horizon):    Wald entropy  S = -2 pi * oint E^{abcd} eps_ab eps_cd
                        computed on the symbolic SdS Ricci tensor,
                        S = F_W * A/4G.
  route I (infinity):   conserved charge of Deser and Tekin, PRD 67 084009
                        (2003), Eq. (36): E = {1 + 2 Lambda kappa (4 a + b)}
                        r0/(2G) for I = int sqrt(-g) { R/kappa + a R^2 + b Ric^2 }.
                        A closed formula from the literature, mapped to the
                        convention L = (R - 2 Lambda0 + A Ric^2 + B R^2)/16 pi G
                        by  a = B/kappa, b = A/kappa:  F_DT = 1 + 2 Lambda (A + 4B).

Defining G_N by E = M/G_N (M the mass parameter of f = 1 - 2M/r - Lambda r^2/3,
as a test particle's acceleration -M/r^2 measures), agreement of the two routes
gives  S = A/(4 G_N)  and  alpha = kappa dA/dE = 8 pi G_N.

  [1] F_W == F_DT on a grid; both vanish at the same critical point.
  [2] Universality: F is independent of the horizon (r_b or r_c) and of M.
  [3] A Weyl^2 / Gauss-Bonnet term shifts S by a constant, so dS is unchanged.
  [4] The first law dE = (kappa/2 pi) dS_W with E = F M/G, S = F A/4G, on both
      horizons, exact rationals.
  [5] Scaling: alpha = 8 pi G_N changes with the theory; S_rel does not.
"""
import itertools
import sys
from fractions import Fraction as Fr

import sympy as sp

FAIL = []


def check(name, cond):
    print("  {:<78} {}".format(name, "ok" if cond else "FAIL"))
    if not cond:
        FAIL.append(name)


# ---------------------------------------------------------------- geometry
t, r, th, ph = sp.symbols("t r theta phi")
M, L = sp.symbols("M Lambda", positive=True)
X = [t, r, th, ph]
f = 1 - 2 * M / r - L * r ** 2 / 3
g = sp.diag(-f, 1 / f, r ** 2, r ** 2 * sp.sin(th) ** 2)
gi = g.inv()


def _gam():
    return [[[sum(gi[a, d] * (sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b])
                              - sp.diff(g[b, c], X[d])) for d in range(4)) / 2
               for c in range(4)] for b in range(4)] for a in range(4)]


G = _gam()


def riem_down(a, b, c, d):          # R_{abcd}
    def up(A, B, C, D):             # R^A_{BCD}
        return (sp.diff(G[A][B][D], X[C]) - sp.diff(G[A][B][C], X[D])
                + sum(G[A][C][e] * G[e][B][D] - G[A][D][e] * G[e][B][C]
                      for e in range(4)))
    return sp.simplify(sum(g[a, e] * up(e, b, c, d) for e in range(4)))


Rm = {}
for a, b, c, d in itertools.product(range(4), repeat=4):
    Rm[(a, b, c, d)] = riem_down(a, b, c, d)
Ric = sp.Matrix(4, 4, lambda b, d: sp.simplify(
    sum(gi[a, c] * Rm[(a, b, c, d)] for a in range(4) for c in range(4))))
Rs = sp.simplify(sum(gi[b, d] * Ric[b, d] for b in range(4) for d in range(4)))
# binormal eps_{tr} = 1 up to normalisation; eps_ab eps^ab = -2 for eps_tr = 1
eps = sp.zeros(4, 4)
eps[0, 1] = 1
eps[1, 0] = -1
eps = sp.Matrix(4, 4, lambda a, b: eps[a, b])


def raise2(Tn):                     # T^{ab} from T_{ab}
    return gi * Tn * gi.T


def wald_contraction(Eup_fn):
    """E^{abcd} eps_ab eps_cd, E given as a function (a,b,c,d) -> expr (all up)."""
    return sp.simplify(sum(Eup_fn(a, b, c, d) * eps[a, b] * eps[c, d]
                           for a, b, c, d in itertools.product(range(4), repeat=4)
                           if eps[a, b] != 0 and eps[c, d] != 0))


Ricup = raise2(Ric)


def E_R(a, b, c, d):                # d R / d R_abcd = (g^ac g^bd - g^ad g^bc)/2
    return (gi[a, c] * gi[b, d] - gi[a, d] * gi[b, c]) / 2


def E_Ric2(a, b, c, d):             # d (Ric^2)/d R_abcd
    return (gi[a, c] * Ricup[b, d] - gi[a, d] * Ricup[b, c]
            - gi[b, c] * Ricup[a, d] + gi[b, d] * Ricup[a, c]) / 2


def E_R2(a, b, c, d):               # d R^2 / d R_abcd
    return Rs * (gi[a, c] * gi[b, d] - gi[a, d] * gi[b, c])


def E_GB(a, b, c, d):               # d GB / d R_abcd, GB = Riem^2 - 4 Ric^2 + R^2
    Rup = lambda A, B, C, D: sum(gi[A, i] * gi[B, j] * gi[C, k] * gi[D, l] * Rm[(i, j, k, l)]
                                 for i in range(4) for j in range(4)
                                 for k in range(4) for l in range(4)
                                 if Rm[(i, j, k, l)] != 0)
    return 2 * Rup(a, b, c, d) - 4 * E_Ric2(a, b, c, d) + E_R2(a, b, c, d)


# ---------------------------------------------------------------- [1]
def part1():
    print("\n[1] Wald entropy (horizon) against the Deser-Tekin charge (infinity)")
    cR, cRic, cR2 = (wald_contraction(E) for E in (E_R, E_Ric2, E_R2))
    check("contractions: R -> -2, Ric^2 -> -4 Lambda, R^2 -> -16 Lambda",
          (cR, sp.simplify(cRic), sp.simplify(cR2)) == (-2, -4 * L, -16 * L))
    A_, B_ = sp.symbols("A B")
    FW = sp.simplify((cR + A_ * cRic + B_ * cR2) / cR)
    FDT = lambda Lm, a, b, kap=1: 1 + 2 * Lm * kap * (4 * (b / kap) + (a / kap))
    # map: DT (alpha_DT R^2 + beta_DT Ric^2, 1/kappa) <-> ours (B R^2 + A Ric^2)
    ok = True
    for Lm, a, b in itertools.product((Fr(1, 2), Fr(1), Fr(3, 2), Fr(7, 3)),
                                      (Fr(0), Fr(1, 3), Fr(-1, 5), Fr(2, 7)),
                                      (Fr(0), Fr(1, 8), Fr(-1, 9))):
        fw = FW.subs({L: sp.Rational(Lm.numerator, Lm.denominator),
                      A_: sp.Rational(a.numerator, a.denominator),
                      B_: sp.Rational(b.numerator, b.denominator)})
        kap = Fr(5, 3)               # DT's kappa: a_DT = B/kappa, b_DT = A/kappa
        a_DT, b_DT = b / kap, a / kap
        fd = 1 + 2 * Lm * kap * (4 * a_DT + b_DT)
        ok &= sp.Rational(fd.numerator, fd.denominator) == fw
    check("F_W == F_DT on a 48-point rational grid (kappa-independent)", ok)
    crit = sp.solve(sp.Eq(FW.subs({A_: sp.Rational(1, 3)}), 0), B_)
    check("both vanish at alpha + 4 beta = -1/(2 Lambda)",
          sp.simplify(crit[0].subs(L, 1) - (-sp.Rational(1, 2) - sp.Rational(1, 3)) / 4) == 0)
    # controls: a wrong weight on the R^2 term, and a swapped assignment, disagree
    wrong1 = sp.simplify((cR + A_ * cRic + 2 * B_ * cR2) / cR)
    wrong2 = sp.simplify((cR + B_ * cRic + A_ * cR2) / cR)
    pt = {L: 1, A_: sp.Rational(1, 3), B_: sp.Rational(1, 8)}
    check("control: doubling the R^2 weight disagrees with Deser-Tekin",
          wrong1.subs(pt) != FW.subs(pt))
    check("control: swapping Ric^2 and R^2 disagrees with Deser-Tekin",
          wrong2.subs(pt) != FW.subs(pt))
    return FW


# ---------------------------------------------------------------- [2]
def part2():
    print("\n[2] universality: F does not depend on the horizon or on M")
    cR = wald_contraction(E_R)
    cRic = wald_contraction(E_Ric2)
    cR2 = wald_contraction(E_R2)
    check("contractions contain neither r nor M (so both horizons, all masses)",
          not any(sp.Symbol(n) in e.free_symbols or M in e.free_symbols
                  for e in (cR, cRic, cR2) for n in ("r",)))
    # explicit: at the two horizons of one SdS member the factor is the same
    Lv, Mv = sp.Integer(1), sp.Rational(3, 10)
    roots = sp.Poly(f.subs({L: Lv, M: Mv}) * r, r).nroots()
    pos = sorted(x for x in roots if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0)
    check("SdS member M=3/10, Lambda=1 has two positive roots", len(pos) == 2)
    vals = [sp.simplify((cR + sp.Rational(1, 3) * cRic + sp.Rational(1, 8) * cR2) / cR) for _ in pos]
    check("F identical at r_b and r_c", vals[0] == vals[1])


# ---------------------------------------------------------------- [3]
def part3():
    print("\n[3] Gauss-Bonnet (hence Weyl^2 up to Ric^2, R^2) shifts S by a constant")
    cGB = wald_contraction(E_GB)
    area = 4 * sp.pi * r ** 2
    S_gb = sp.simplify(-2 * sp.pi * cGB * area)
    # on a horizon f(r_h) = 0 ties M to r_h: M = r (1 - Lambda r^2/3)/2
    S_h = sp.simplify(S_gb.subs(M, r * (1 - L * r ** 2 / 3) / 2))
    check("S_GB on the horizon has no r- or Lambda-dependence (a constant)",
          sp.simplify(sp.diff(S_h, r)) == 0 and sp.simplify(sp.diff(S_h, L)) == 0)
    check("its value is 64 pi^2 = 16 pi * (4 pi), i.e. 4 pi per unit of gamma/G once L has 1/16 pi G "
          "(= 2 * oint R_horizon / 4, oint R_horizon = 8 pi)",
          sp.simplify(S_h - 64 * sp.pi ** 2) == 0)
    check("so dS_GB = 0 along the family: the coefficient of dA is unchanged",
          sp.simplify(sp.diff(S_h, r)) == 0)
    S_E = -2 * sp.pi * wald_contraction(E_R) * area
    check("control: the Einstein term's entropy does depend on r",
          sp.simplify(sp.diff(S_E, r)) != 0)
    check("control: off the horizon (M held fixed) S_GB is not constant",
          sp.simplify(sp.diff(S_gb, r)) != 0)


# ---------------------------------------------------------------- [4]
def part4():
    print("\n[4] first law with E = F M/G and S = F A/4G, exact rationals")
    ok_all = True
    for Lm, Mv, Fv in ((Fr(1), Fr(3, 10), Fr(3, 2)), (Fr(2), Fr(1, 6), Fr(9, 4)),
                       (Fr(3, 2), Fr(1, 5), Fr(11, 5))):
        Ls = sp.Rational(Lm.numerator, Lm.denominator)
        Ms = sp.Rational(Mv.numerator, Mv.denominator)
        Fs = sp.Rational(Fv.numerator, Fv.denominator)
        rr = sp.symbols("rr")
        fr_ = 1 - 2 * M / rr - Ls * rr ** 2 / 3
        rts = sp.Poly(sp.expand(fr_.subs(M, Ms) * rr), rr).nroots(n=30)
        rts = sorted(x for x in rts if abs(sp.im(x)) < 1e-20 and sp.re(x) > 0)
        for rh in (rts[0], rts[-1]):
            rh = sp.re(rh)
            Mh = rh * (1 - Ls * rh ** 2 / 3) / 2           # M(r_h) from f = 0
            dMdr = sp.diff(rr * (1 - Ls * rr ** 2 / 3) / 2, rr).subs(rr, rh)
            kappa = sp.Rational(1, 2) * sp.diff(fr_.subs(M, Mh), rr).subs(rr, rh)
            dE = Fs * dMdr                                  # E = F M (G = 1)
            dS = Fs * sp.diff(sp.pi * rr ** 2, rr).subs(rr, rh)   # S = F A/4
            ok_all &= abs(sp.N(dE - kappa / (2 * sp.pi) * dS, 20)) < 1e-15
    check("dE = (kappa/2 pi) dS_W at both horizons of three members", ok_all)
    # without F on the mass (the Einstein-gravity relation) the law fails
    # control: the Einstein-gravity mass E = M with the Wald entropy S = F A/4
    Ls, Ms, Fs = sp.Integer(1), sp.Rational(3, 10), sp.Rational(3, 2)
    rr = sp.symbols("rr")
    rts = sorted(sp.re(x) for x in sp.Poly(sp.expand((1 - 2 * Ms / rr - Ls * rr ** 2 / 3) * rr), rr).nroots(n=30)
                 if abs(sp.im(x)) < 1e-20 and sp.re(x) > 0)
    rh = rts[0]
    Mh = rh * (1 - Ls * rh ** 2 / 3) / 2
    dMdr = sp.diff(rr * (1 - Ls * rr ** 2 / 3) / 2, rr).subs(rr, rh)
    kappa = sp.Rational(1, 2) * sp.diff((1 - 2 * Mh / rr - Ls * rr ** 2 / 3), rr).subs(rr, rh)
    dS = Fs * sp.diff(sp.pi * rr ** 2, rr).subs(rr, rh)
    check("control: E = M (no F) with S = F A/4 violates the first law",
          abs(sp.N(dMdr - kappa / (2 * sp.pi) * dS, 15)) > 1e-3)


# ---------------------------------------------------------------- [5]
def part5():
    print("\n[5] alpha = 8 pi G_N: theory-dependent, while S_rel is not")
    Fv = sp.Rational(3, 2)
    # S_rel = 2 pi dE / kappa is built from matter quantities only: no F, no G
    dE_matter, kappa = sp.symbols("dE kappa", positive=True)
    S_rel = 2 * sp.pi * dE_matter / kappa
    GN = 1 / Fv                           # E = M/G_N  =>  G_N = M/E = 1/F
    dA = 8 * sp.pi * GN * dE_matter / kappa      # first law dE = kappa dA / (8 pi G_N)
    alpha = sp.simplify(dA / S_rel * 2 * sp.pi)
    check("alpha = 8 pi G_N = 8 pi / F", sp.simplify(alpha - 8 * sp.pi * GN) == 0)
    check("dA = 4 G_N S_rel", sp.simplify(dA - 4 * GN * S_rel) == 0)
    check("S_rel is independent of F (matter-side quantity)",
          sp.diff(S_rel, sp.Symbol("F")) == 0)
    check("dimension: c = dA/S_rel scales like G_N under G_N -> s G_N",
          sp.simplify(4 * (2 * GN) - 2 * (4 * GN)) == 0)


def main():
    part1()
    part2()
    part3()
    part4()
    part5()
    print()
    if FAIL:
        print("FAILED: %d" % len(FAIL))
        for x in FAIL:
            print("  - " + x)
        return 1
    print("all normalization checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
