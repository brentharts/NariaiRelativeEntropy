# Nariai Relative Entropy

Companion repository for the paper

> **Relative Entropy on the Nariai Horizon: A Finite-Area Completion of the
> Entropic Derivation of the Semiclassical Einstein Equations**
> (B. S. Hartshorn, 2026, draft)

which extends Dorau & Much, *From Quantum Relative Entropy to the
Semiclassical Einstein Equations*, Phys. Rev. Lett. **136**, 091602 (2026)
[arXiv:2510.24491], from the local Rindler wedge to the Nariai spacetime
dS2 x S2 -- the degenerate limit of Schwarzschild-de Sitter in which the
black-hole and cosmological horizons merge.

## Why Nariai

* **Finiteness.** The bifurcation surface is a compact two-sphere of area
  4*pi/Lambda, so the relative entropy, the area, and the identification
  S_rel = dA/4 all involve finite quantities. The infinite-area caveat of the
  Rindler construction disappears, and the ledger equation
  `S = pi/Lambda + S_rel` relates finite numbers.
* **Unambiguous temperature.** The Killing normalization is fixed by the
  geometry (Bousso-Hawking), removing the Rindler boost-rescaling ambiguity.
* **Dynamics.** Positivity of relative entropy becomes a selection rule on the
  branches of the Nariai (anti-)evaporation instability: within the coherent
  sector only the pierced horizon cut can grow; anti-evaporation is a
  certificate of non-coherent horizon flux.

## Key exact result (Appendix B)

With `U = -exp(-kappa*u)/kappa`, the modular weight combines with the boost
Jacobian to leave the boost measure invariant:

```
(-U) |d_U chi|^2 dU  =  (1/kappa) |d_u chi|^2 du
```

Note the derivative factors are needed. The weight and the measure *alone*
do not cancel: `(-U)(dU/du) = exp(-2*kappa*u)/kappa`, not `1/kappa`. Since
`d_U = exp(kappa*u) d_u`, a squared first derivative supplies the missing
`exp(+2*kappa*u)`. In exponents of `exp(kappa*u)` the three factors are
`-1, +2, -1`, and the identity is that they sum to zero
(`NariaiFacts.lean: modular_weight_balances`). Hence the interference
functional

```
I2 = Int (-U) (d_U chi)^2 dU dvol
```

is the zero-total-boost-frequency Fourier component of `(d_u chi)^2` and
**vanishes identically** for boost-positive-frequency packets. Consequences:

1. The squeezed-state relative entropy takes the thermal first-law form
   `S_rel = beta * dE_boost`, `beta = 2*pi/sqrt(Lambda)`.
2. Local flux-negativity windows (the engine of anti-evaporation) integrate
   to zero against the modular weight: negativity is a zero-sum reallocation.
3. `|I2|/I1` is a wedge-locality detector: 0 for boost-adapted data, 1 for
   boost-blind data. Between those endpoints it is *graded*, not binary. For
   a negative-frequency admixture of weight `eps` the closed form is

   ```
   |I2|/I1 = 2*eps/(1 + eps^2) = tanh(2*artanh(eps))
   ```

   so the squeeze-shape inequality `|I2|/I1 <= tanh(r)` reads
   `r >= 2*artanh(eps)`: a minimum squeeze, diverging as `eps -> 1`. That
   divergence is the sharp form of "no wedge-local squeeze exists for real
   data". Verified to 1e-15 at six values of `eps`; the two integer facts
   underneath (the ratio is bounded by one, and saturates only at `eps = 1`)
   are proved in `NariaiFacts.lean`.

## Contents

| File | Purpose |
|---|---|
| `supplementary.tex` | Supplementary material: TikZ diagrams + numerical figures |
| `appendix_b_verification.py` | Numerical checks for Appendix B (I2 vanishing, negativity windows, one-mode sum rule) |
| `nariailean.py` | Emits `NariaiFacts.lean`: the finite algebra behind Appendix B |
| `NariaiFacts.lean` | Mathlib-free Lean 4, kernel-checked, no admitted proofs |
| `generate_figures.py` | Generates Figures S1-S3 (matplotlib) into `figures/` |
| `Makefile` | Build system |

## Quickstart

```sh
make verify        # run the numerical checks (numpy)
make selftest      # pass/fail assertions on every claim
make convergence   # grid and truncation study for the I2 vanishing
make lean          # regenerate and kernel-check NariaiFacts.lean
make figures       # regenerate Figures S1-S3 (matplotlib)
make               # build supplementary.pdf (needs pdflatex + TikZ)
```

Requirements: Python >= 3.10 with `numpy` and `matplotlib`; a TeX
distribution with TikZ (TeX Live recommended).

## Verified numbers

From `make verify` (kappa = 1):

* `|I2|/I1 = 1.3e-15` for a boost-positive-frequency packet (prediction: 0).
* Boost/affine cross-check: `1.6e-7`, falling as `2.6e-6 -> 6.4e-7 -> 1.6e-7`
  under grid refinement. This differentiates `chi` with respect to `U`
  numerically on the `U` grid, supplying no analytic Jacobian, so it is a
  check that can fail -- and a deliberately wrong `kappa` in the coordinate
  map does break it. (An earlier version reported `2e-16` here, but that
  computation reduced symbolically to the boost-coordinate integrand and
  could not fail; the number was floating-point noise.)
* `|I2|/I1 = 1.000000` for a real (boost-blind) packet, computed through the
  same complex code path as the boost-adapted case rather than asserted.
* Detector closed form `2*eps/(1+eps^2)` reproduced to `1e-15` at six values
  of `eps`.
* Squeezed flux at `(r, theta) = (0.5, pi)`: negative on 51% of the packet
  support, funded depth 21% of the positive part, total positive. The total
  is `theta`-independent to `6e-15`, which is the zero-sum theorem.
* One-mode sum rule `S_rel = beta * dE` verified to <= 2e-15 at three
  parameter points, matching `beta * w * sinh(r)^2 * coth(beta*w/2)`. The
  reference state returns `S_th = 0` to `3e-16`, and the trace gap matches
  the Lean form `nu*(s-t)^2` to `2e-16`.

### How generic is the `1e-15`?

Not very, and `make convergence` says so. The vanishing is exact only for
spectra bounded away from zero on the whole line, and the numerics have two
finite-size errors:

| effect | truncated Gaussian | smooth-at-zero spectrum |
|---|---|---|
| carrier `w0 = 3.0` | `1.3e-15` | `5.9e-17` |
| carrier `w0 = 1.5` | `3.6e-08` | `2.0e-14` |
| carrier `w0 = 1.0` | `1.4e-07` | `1.0e-12` |

The degradation is the hard truncation of the Gaussian at `w_min > 0`, not
the theorem: a spectrum vanishing smoothly at the origin loses six orders
less. The residual is the finite `u` window, and it converges cleanly --
`7.5e-08 -> 5.7e-12 -> 2.0e-14` as the window widens from `+-30` to `+-120`.

## Machine-checked fragments

`NariaiFacts.lean` (17 theorems, Mathlib-free, kernel-accepted, no `sorryAx`
and no `Classical.choice`) proves the parts of Appendix B that are finite
algebra rather than analysis:

* the modular-weight exponent balance, and the failure of the naive form;
* the mechanism of `I2 = 0` as a statement about a **sumset** -- `I2` is the
  zero-frequency component of `(d_u chi)^2`, whose spectrum is the sumset of
  the packet's, so `I2 = 0` exactly when no two frequencies sum to zero.
  Proved for arbitrary finite spectra, both the positive case (never) and
  the symmetric case (always);
* the one-mode sum rule as ring identities: the determinant is unchanged
  (hence `dS_vN = 0`), and the trace gap factors as `nu*(s-t)^2` (hence the
  closed form, and non-negativity of the energy cost);
* the detector bound and its saturation condition.

The analysis stays in the paper. What Lean settles is the algebra.

## Status and open problems

* (i) Squeeze-shape inequality: **closed** (exact vanishing theorem).
* (ii) dS = 0 for the quadratic generator: verified at the one-mode level;
  the type-III domain-theoretic lift remains open.
* (iii) Backreaction with sign-indefinite flux (teleological horizon):
  the open dynamical frontier, constrained by the zero-sum theorem.

## Papers
- https://doi.org/10.5281/zenodo.21303604
- https://doi.org/10.5281/zenodo.21382744

