#!/usr/bin/env python3
r"""Numerical verification for Appendix B of
"Relative Entropy on the Nariai Horizon".

Three checks, kappa = 1 throughout (Nariai: kappa = sqrt(Lambda)):

  (B1) For boost-positive-frequency horizon packets the interference
       functional I2 vanishes in the continuum, so the squeeze-shape
       inequality |I2|/I1 <= tanh r holds trivially.
  (B2) The squeezed flux <T_UU> has genuinely negative local windows while
       the modular-weighted total stays positive.
  (B3) One-mode Gaussian check of the sum rule S_rel = beta * Delta E.

What changed, and why
---------------------
The previous version of this file reported two numbers that could not fail,
and both were quoted in the README as evidence.

* The "affine U cross-check" defined `dchi_dU = dchi_du * exp(kappa*u)` and
  then integrated `(-U) |dchi_dU|^2 exp(-kappa*u) du`.  Substituting
  `(-U) = exp(-kappa*u)/kappa` the exponentials cancel exactly and the
  integrand reduces to `|dchi_du|^2 / kappa` -- symbolically identical to
  I1's.  The reported agreement of 2e-16 was floating-point noise from
  summing the same numbers in a different order, and the Jacobian it claimed
  to verify had been imposed by hand one line earlier.  It is replaced here
  by an INDEPENDENT check: chi is differentiated with respect to U
  numerically, on the U grid, with no analytic Jacobian supplied, and the
  functional is integrated over U.  That agrees to ~1e-7 and improves under
  refinement, which is what a real check looks like.

* The (B1') contrast case assigned `I1R` and `I2R` the same expression and
  printed their ratio as `1.000000`.  It is now computed through the same
  complex code path as (B1), so the equality is measured; and a family
  interpolating between boost-adapted and boost-blind data is swept, which
  tests the "binary detector" claim instead of restating its endpoints.

* A convergence study was added.  The headline 1e-15 is not generic: it
  degrades to ~1e-7 when the packet's spectral support approaches zero.  Two
  distinct finite-size errors are separated below -- the hard truncation of
  the spectrum at w_min, and the finite u window -- and neither is a failure
  of the theorem, which is exact for spectra bounded away from zero on the
  whole line.  Reporting only the best case would have been the same mistake
  as the two identities.

Also: the README states the modular-weight identity as `(-U) dU = du/kappa`.
Taken literally that is false -- `(-U)(dU/du) = exp(-2 kappa u)/kappa`.  The
true identity needs the derivative factors that accompany it in the
integrand.  `check_modular_weight_exponents` below states the correct form,
and NariaiFacts.lean proves the exponent balance.

    python3 appendix_b_verification.py
    python3 appendix_b_verification.py --selftest
    python3 appendix_b_verification.py --convergence
"""

import sys

import numpy as np

KAPPA = 1.0

# Exponents of e^{kappa u} carried by the three factors of the integrand.
# Mirrors nariailean.WEIGHT_EXP / DERIV_EXP / JACOBIAN_EXP; the Lean file
# proves they sum to zero and that dropping the derivative factor leaves -2.
WEIGHT_EXP, DERIV_EXP, JACOBIAN_EXP = -1, 2, -1


# ----------------------------------------------------------------------
# packets
# ----------------------------------------------------------------------

def gaussian_spectrum(w0=3.0, sg=0.7, wmin=1e-4, wmax=12.0, nw=2400):
    """A Gaussian in frequency, hard-truncated at wmin > 0.

    The truncation is the whole reason this is not exactly positive
    frequency: the Gaussian tail is cut, and the cut is an edge.  For a
    well-separated packet the edge is at e^-20 and invisible; for a broad or
    low-carrier one it is not, which the convergence study shows.
    """
    w = np.linspace(wmin, wmax, nw)
    return w, np.exp(-(w - w0) ** 2 / (2 * sg ** 2))


def smooth_spectrum(w0=3.0, wmax=20.0, nw=2400, power=4):
    """A spectrum vanishing to high order at w = 0, with no edge.

    `(w/w0)^p exp(-4w/w0)` is zero at the origin along with its first p
    derivatives, so truncating the grid at w = 0 introduces no
    discontinuity.  Comparing this against the truncated Gaussian is what
    separates the spectral-edge error from the u-window error.
    """
    w = np.linspace(0.0, wmax, nw)
    return w, (w / w0) ** power * np.exp(-4.0 * w / w0)


def build_packet(w, aw, u, chunk=4000):
    """chi(u) and d_u chi(u) for chi = int dw a(w) e^{-i w u}."""
    dw = w[1] - w[0]
    c0, c1 = aw * dw, (-1j * w * aw) * dw
    chi = np.empty(u.shape, dtype=complex)
    dchi = np.empty(u.shape, dtype=complex)
    for i0 in range(0, len(u), chunk):
        blk = u[i0:i0 + chunk]
        phase = np.exp(-1j * np.outer(blk, w))
        chi[i0:i0 + chunk] = phase @ c0
        dchi[i0:i0 + chunk] = phase @ c1
    return chi, dchi


# ----------------------------------------------------------------------
# (B1) the functionals
# ----------------------------------------------------------------------

def functionals_boost(dchi, du, kappa=KAPPA):
    """I1 and I2 in boost coordinates."""
    I1 = np.trapezoid(np.abs(dchi) ** 2, dx=du) / kappa
    I2 = np.trapezoid(dchi ** 2, dx=du) / kappa
    return I1, I2


def functionals_affine(chi, u, kappa=KAPPA, trim=1.0):
    """I1 and I2 computed in affine U coordinates, INDEPENDENTLY.

    chi is differentiated with respect to U by `np.gradient(chi, U)` on the
    (highly non-uniform) U grid, and the functional is integrated over U.
    No analytic Jacobian is used anywhere, so this can disagree with the
    boost-coordinate answer -- and if the change of variables were wrong, it
    would.  The endpoints are trimmed because a one-sided difference there
    is first-order accurate.
    """
    U = -np.exp(-kappa * u) / kappa
    dchi_dU = np.gradient(chi, U)
    m = (u > u[0] + trim) & (u < u[-1] - trim)
    I1 = np.trapezoid((-U[m]) * np.abs(dchi_dU[m]) ** 2, U[m])
    I2 = np.trapezoid((-U[m]) * dchi_dU[m] ** 2, U[m])
    return I1, I2


def affine_crosscheck(window=10.0, nu_pts=40001, kappa=KAPPA,
                      kappa_map=None, **kw):
    """Compare I1 in boost and affine coordinates on a COMMON grid.

    This needs its own, narrower window.  U = -e^{-ku}/k stretches the grid
    exponentially, so at u = -40 the spacing is ~1e17 times the spacing at
    u = 0 and a finite difference there is meaningless -- the check would
    fail for a reason that has nothing to do with the change of variables.
    The packet has support within a few units of the origin, so a window of
    +-10 contains it and keeps the U grid usable.  Both integrals are taken
    over the same trimmed interior, so the comparison is like for like.
    """
    w, aw = gaussian_spectrum(**kw)
    u = np.linspace(-window, window, nu_pts)
    chi, dchi = build_packet(w, aw, u)
    m = (u > u[0] + 1.0) & (u < u[-1] - 1.0)
    I1_boost = np.trapezoid(np.abs(dchi[m]) ** 2, dx=u[1] - u[0]) / kappa
    # kappa_map corrupts ONLY the coordinate map, leaving the boost-side
    # integral alone.  Passing a wrong kappa to both would rescale the two
    # sides together and they would still agree -- a negative control that
    # controls nothing, which is the failure mode this whole file is about.
    I1_affine, I2_affine = functionals_affine(chi, u, kappa_map or kappa)
    return I1_boost, I1_affine, abs(I1_boost - I1_affine) / I1_boost


def check_modular_weight_exponents():
    """The identity the integrand actually satisfies, in exponents.

    (-U) ~ e^{-ku}, |d_U chi|^2 ~ e^{+2ku}, dU ~ e^{-ku}.  They cancel.
    Dropping the derivative factor leaves e^{-2ku}, which is the README's
    literal claim and is wrong.
    """
    return (WEIGHT_EXP + DERIV_EXP + JACOBIAN_EXP,
            WEIGHT_EXP + JACOBIAN_EXP)


def ratio_for(spectrum=gaussian_spectrum, window=40.0, nu_pts=32001,
              **kw):
    """|I2|/I1 for one packet and one grid."""
    w, aw = spectrum(**kw)
    u = np.linspace(-window, window, nu_pts)
    _, dchi = build_packet(w, aw, u)
    I1, I2 = functionals_boost(dchi, u[1] - u[0])
    return abs(I2) / I1


# ----------------------------------------------------------------------
# (B1') the detector, swept rather than asserted
# ----------------------------------------------------------------------

def mixed_packet_ratio(eps, w0=3.0, sg=0.7, window=40.0, nu_pts=32001):
    """|I2|/I1 for a packet with a negative-frequency admixture of weight eps.

    eps = 0 is boost-adapted (positive spectrum only); eps = 1 is a real
    packet, whose spectrum is symmetric.  The README calls |I2|/I1 a binary
    detector, 0 or 1.  The endpoints are binary; this measures what happens
    between them, which the previous version could not, because it asserted
    both endpoints instead of computing either.
    """
    w, aw = gaussian_spectrum(w0=w0, sg=sg)
    u = np.linspace(-window, window, nu_pts)
    _, dpos = build_packet(w, aw, u)
    dneg = np.conj(dpos)                 # the negative-frequency mirror
    dchi = dpos + eps * dneg
    I1, I2 = functionals_boost(dchi, u[1] - u[0])
    return abs(I2) / I1


def detector_closed_form(eps):
    """|I2|/I1 = 2 eps / (1 + eps^2) for an admixture of weight eps.

    Derived, not fitted.  With d_u chi = P + eps conj P and int P^2 = 0 from
    the vanishing theorem, I1 = (1 + eps^2) a and I2 = 2 eps a where
    a = int |P|^2, so the ratio is this and the packet normalisation drops
    out entirely.  Equivalently 2 eps/(1 + eps^2) = tanh(2 artanh eps), so
    the squeeze-shape inequality |I2|/I1 <= tanh r becomes r >= 2 artanh eps:
    an admixture of weight eps demands at least that much squeeze, and the
    demand diverges as eps -> 1, which is the sharp form of "no wedge-local
    squeeze exists for real data".

    NariaiFacts.lean proves the two integer facts underneath: the ratio is
    bounded by one, and it saturates exactly when eps = 1.
    """
    return 2.0 * eps / (1.0 + eps ** 2)


def min_squeeze_for(eps):
    """r_min = 2 artanh eps, the least squeeze representing that admixture."""
    return 2.0 * np.arctanh(eps) if eps < 1.0 else np.inf


# ----------------------------------------------------------------------
# (B2) negativity windows
# ----------------------------------------------------------------------

def flux_profile(dchi, r, theta, kappa=KAPPA):
    """Modular-weighted flux density in boost time."""
    return (2 * np.sinh(r) ** 2 * np.abs(dchi) ** 2
            - np.sinh(2 * r) * np.real(np.exp(-1j * theta) * dchi ** 2)) / kappa


def flux_summary(dchi, du, r, theta):
    t = flux_profile(dchi, r, theta)
    total = np.trapezoid(t, dx=du)
    neg = -np.trapezoid(np.where(t < 0, t, 0.0), dx=du)
    pos = np.trapezoid(np.where(t > 0, t, 0.0), dx=du)
    return dict(total=total, depth=neg, pos=pos, frac=float(np.mean(t < 0)))


# ----------------------------------------------------------------------
# (B3) the one-mode sum rule
# ----------------------------------------------------------------------

def svn(nu):
    a, b = (nu + 1) / 2, (nu - 1) / 2
    return a * np.log(a) - (b * np.log(b) if b > 0 else 0.0)


def one_mode(w_m, r, kappa=KAPPA):
    beta = 2 * np.pi / kappa
    nu = 1.0 / np.tanh(beta * w_m / 2)
    sig_th = nu * np.eye(2)
    S = np.diag([np.exp(r), np.exp(-r)])
    sig_sq = S @ sig_th @ S.T
    nu_sq = np.sqrt(np.linalg.det(sig_sq))
    E_th = w_m * np.trace(sig_th) / 4
    E_sq = w_m * np.trace(sig_sq) / 4
    lnZ = -np.log(2 * np.sinh(beta * w_m / 2))
    S_rel = beta * E_sq + lnZ - svn(nu_sq)
    return dict(beta=beta, nu=nu, nu_sq=nu_sq,
                dS_vN=svn(nu_sq) - svn(nu),
                S_rel=S_rel, beta_dE=beta * (E_sq - E_th),
                # the reference state must give exactly zero; the previous
                # version computed this and never printed it
                S_th=beta * E_th + lnZ - svn(nu),
                analytic=beta * w_m * nu * np.sinh(r) ** 2,
                # the Lean identity: trace gap = nu (s - t)^2 with s = e^r
                trace_gap=float(np.trace(sig_sq) - np.trace(sig_th)),
                lean_form=nu * (np.exp(r) - np.exp(-r)) ** 2)


# ----------------------------------------------------------------------
# reporting
# ----------------------------------------------------------------------

def report():
    w, aw = gaussian_spectrum()
    u = np.linspace(-40.0, 40.0, 32001)
    du = u[1] - u[0]
    chi, dchi = build_packet(w, aw, u)

    bal, naive = check_modular_weight_exponents()
    print('== (B0) the modular weight, in exponents of e^{kappa u} ==')
    print(' (-U) %+d,  |d_U chi|^2 %+d,  dU %+d   -> sum = %d  (must be 0)'
          % (WEIGHT_EXP, DERIV_EXP, JACOBIAN_EXP, bal))
    print(' weight and measure alone   -> sum = %d, so "(-U) dU = du/kappa"'
          % naive)
    print(' is false as written; the derivative factor is what cancels it.')

    I1, I2 = functionals_boost(dchi, du)
    print('\n== (B1) boost-positive-frequency packet ==')
    print(' I1 (boost coords)            = %.6e' % I1)
    print(' |I2|/I1                      = %.3e   (continuum: 0)'
          % (abs(I2) / I1))
    print(' affine U cross-check, on a common +-10 window:')
    for n in (20001, 40001, 80001):
        b, a, rel = affine_crosscheck(nu_pts=n)
        print('   n=%-8d I1_boost = %.6e  I1_affine = %.6e  rel = %.2e'
              % (n, b, a, rel))
    print('   (numerical d/dU on the U grid; no analytic Jacobian used, and')
    print('    the error falls under refinement, so it is a real check.)')

    print('\n== (B1\') the detector, swept ==')
    print(' %-8s %-14s %-14s %-10s %s'
          % ('eps', '|I2|/I1', '2e/(1+e^2)', 'diff', 'r_min = 2 artanh e'))
    for eps in (0.0, 0.05, 0.2, 0.5, 0.8, 1.0):
        m, c = mixed_packet_ratio(eps), detector_closed_form(eps)
        print(' %-8.2f %-14.8f %-14.8f %-10.2e %.4f'
              % (eps, m, c, abs(m - c), min_squeeze_for(eps)))
    print(' The endpoints are 0 and 1, so the README\'s binary reading is')
    print(' right about them -- but between them the detector is graded, and')
    print(' the grading is exactly tanh(2 artanh eps).  So the content is a')
    print(' minimum squeeze r >= 2 artanh eps, diverging as eps -> 1.')

    print('\n== (B2) squeezed flux: local negativity, global positivity ==')
    rows = []
    for r, th in [(0.5, 0.0), (0.5, np.pi), (1.0, np.pi / 2)]:
        s = flux_summary(dchi, du, r, th)
        rows.append((r, th, s))
        print(' r=%.1f th=%5.2f:  S_rel/2pi = % .5e  D(W) = %.5e  '
              'D/pos = %.3f  neg support = %4.1f%%  total>0: %s'
              % (r, th, s['total'], s['depth'], s['depth'] / s['pos'],
                 100 * s['frac'], s['total'] > 0))
    a, b = rows[0][2]['total'], rows[1][2]['total']
    print(' theta-independence of the total: |diff|/total = %.2e' %
          (abs(a - b) / abs(a)))
    print(' (that is the zero-sum theorem: the theta term is the I2 term,')
    print('  and it integrates to zero, so only the support moves.)')

    print('\n== (B3) one-mode sum-rule check ==')
    for w_m, r in [(0.5, 0.3), (1.0, 0.8), (2.0, 1.5)]:
        d = one_mode(w_m, r)
        print(' w=%.1f r=%.1f: nu_sq-nu = %.2e  dS_vN = %.2e  S_rel = %.6f'
              % (w_m, r, d['nu_sq'] - d['nu'], d['dS_vN'], d['S_rel']))
        print('            beta*dE = %.6f   |S_rel-beta*dE| = %.2e'
              % (d['beta_dE'], abs(d['S_rel'] - d['beta_dE'])))
        print('            analytic beta*w*nu*sinh^2 r = %.6f' % d['analytic'])
        print('            reference state S_th = %.2e (must be 0)' % d['S_th'])
        print('            Lean form nu(s-t)^2 vs trace gap: %.2e'
              % abs(d['trace_gap'] - d['lean_form']))


def convergence():
    """Separate the two finite-size errors in the I2 vanishing."""
    print('== convergence of |I2|/I1 ==')
    print('\n grid refinement, well-separated packet (w0=3.0, sg=0.7)')
    for lbl, kw in [('baseline', {}),
                    ('u-grid x2', {'nu_pts': 64001}),
                    ('u-window +-80', {'window': 80.0}),
                    ('w-grid x2', {'nw': 4800})]:
        print('   %-34s %.3e' % (lbl, ratio_for(**kw)))

    print('\n spectral edge: truncated Gaussian vs smooth-at-zero spectrum')
    print('   %-14s %-16s %s' % ('carrier w0', 'truncated', 'smooth'))
    for w0, sg in [(3.0, 0.7), (1.5, 1.0), (1.0, 2.0)]:
        g = ratio_for(w0=w0, sg=sg)
        s = ratio_for(spectrum=smooth_spectrum, w0=w0)
        print('   %-14.1f %-16.3e %.3e' % (w0, g, s))
    print('   The truncated Gaussian degrades by orders of magnitude as its')
    print('   support nears w=0; the smooth spectrum, which has no edge,')
    print('   degrades far less.  So that error is the truncation, not the')
    print('   theorem -- which is exact for spectra bounded away from zero.')

    print('\n u-window, low-carrier smooth packet (the residual error)')
    for W, n in [(30.0, 48001), (60.0, 96001), (120.0, 192001)]:
        print('   window +-%-6.0f n=%-8d %.3e'
              % (W, n, ratio_for(spectrum=smooth_spectrum, w0=0.5,
                                 window=W, nu_pts=n)))
    print('   A broad packet in frequency is narrow nowhere in u, so the')
    print('   finite window clips it.  Widening the window is what fixes it.')


# ----------------------------------------------------------------------

def selftest():
    failures = []

    def check(label, ok):
        print('  %-58s %s' % (label, 'ok' if ok else 'FAIL'))
        if not ok:
            failures.append(label)

    w, aw = gaussian_spectrum()
    u = np.linspace(-40.0, 40.0, 32001)
    du = u[1] - u[0]
    chi, dchi = build_packet(w, aw, u)
    I1, I2 = functionals_boost(dchi, du)

    print('the modular weight')
    bal, naive = check_modular_weight_exponents()
    check('the three exponents cancel', bal == 0)
    check('weight and measure alone do not', naive != 0)
    check('and leave exactly -2', naive == -2)

    print('(B1) the vanishing, and a cross-check that can fail')
    check('I1 is positive', I1 > 0)
    check('|I2|/I1 is below 1e-12', abs(I2) / I1 < 1e-12)
    b1, a1, rel1 = affine_crosscheck(nu_pts=20001)
    b2, a2, rel2 = affine_crosscheck(nu_pts=80001)
    check('the independent affine I1 agrees to 1e-5', rel2 < 1e-5)
    check('and the agreement improves under refinement', rel2 < rel1)
    # a check that cannot fail is worth nothing, so confirm this one is
    # sensitive: the wrong kappa in the change of variables must break it
    bw, aw_, relbad = affine_crosscheck(nu_pts=40001, kappa_map=2.0)
    check('a wrong kappa in the change of variables breaks it',
          relbad > 1e-3)

    print("(B1') the detector")
    check('boost-adapted end is 0', mixed_packet_ratio(0.0) < 1e-12)
    check('real end is 1', abs(mixed_packet_ratio(1.0) - 1.0) < 1e-9)
    mid = [mixed_packet_ratio(e) for e in (0.05, 0.2, 0.5)]
    check('the sweep is monotone between them',
          all(mid[i] < mid[i + 1] for i in range(len(mid) - 1)))
    check('and takes values strictly inside (0,1)',
          all(0 < m < 1 for m in mid))
    check('it matches the closed form 2eps/(1+eps^2) to 1e-12',
          all(abs(mixed_packet_ratio(e) - detector_closed_form(e)) < 1e-12
              for e in (0.0, 0.05, 0.2, 0.5, 0.8, 1.0)))
    # the two facts NariaiFacts.lean proves, checked numerically here too
    check('the detector never exceeds one (Lean: detector_bounded)',
          all(detector_closed_form(e) <= 1.0 + 1e-15
              for e in np.linspace(0, 1, 101)))
    check('and saturates only at eps=1 (Lean: detector_saturates_iff_equal)',
          all(detector_closed_form(e) < 1.0 - 1e-12
              for e in np.linspace(0, 0.99, 100)))
    check('the closed form is tanh(2 artanh eps)',
          all(abs(detector_closed_form(e) - np.tanh(min_squeeze_for(e))) < 1e-12
              for e in (0.05, 0.2, 0.5, 0.8)))

    print('(B2) negativity is local, positivity is global')
    for r, th in [(0.5, 0.0), (0.5, np.pi), (1.0, np.pi / 2)]:
        s = flux_summary(dchi, du, r, th)
        check('r=%.1f th=%.2f: total positive' % (r, th), s['total'] > 0)
        check('r=%.1f th=%.2f: a genuine negative window exists' % (r, th),
              s['depth'] > 0 and s['frac'] > 0.01)
    t0 = flux_summary(dchi, du, 0.5, 0.0)['total']
    tp = flux_summary(dchi, du, 0.5, np.pi)['total']
    check('the total is theta-independent to 1e-9',
          abs(t0 - tp) / abs(t0) < 1e-9)

    print('(B3) the sum rule')
    for w_m, r in [(0.5, 0.3), (1.0, 0.8), (2.0, 1.5)]:
        d = one_mode(w_m, r)
        check('w=%.1f r=%.1f: symplectic eigenvalue invariant' % (w_m, r),
              abs(d['nu_sq'] - d['nu']) < 1e-12)
        check('w=%.1f r=%.1f: dS_vN vanishes' % (w_m, r),
              abs(d['dS_vN']) < 1e-12)
        check('w=%.1f r=%.1f: S_rel = beta dE' % (w_m, r),
              abs(d['S_rel'] - d['beta_dE']) < 1e-12)
        check('w=%.1f r=%.1f: matches the closed form' % (w_m, r),
              abs(d['S_rel'] - d['analytic']) < 1e-9)
        check('w=%.1f r=%.1f: the reference state gives zero' % (w_m, r),
              abs(d['S_th']) < 1e-12)
        check('w=%.1f r=%.1f: trace gap is nu(s-t)^2, as proved in Lean'
              % (w_m, r), abs(d['trace_gap'] - d['lean_form']) < 1e-9)

    print()
    if failures:
        print('%d failure(s): %s' % (len(failures), ', '.join(failures)))
    else:
        print('appendix_b_verification: all checks pass.')
    return len(failures)


if __name__ == '__main__':
    if '--selftest' in sys.argv:
        sys.exit(1 if selftest() else 0)
    elif '--convergence' in sys.argv:
        convergence()
    else:
        report()
