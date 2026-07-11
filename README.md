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

With `U = -exp(-kappa*u)/kappa`, the modular weight is exactly the boost
Jacobian: `(-U) dU = du / kappa`. Hence the interference functional

```
I2 = Int (-U) (d_U chi)^2 dU dvol
```

is the zero-total-boost-frequency Fourier component of `(d_u chi)^2` and
**vanishes identically** for boost-positive-frequency packets. Consequences:

1. The squeezed-state relative entropy takes the thermal first-law form
   `S_rel = beta * dE_boost`, `beta = 2*pi/sqrt(Lambda)`.
2. Local flux-negativity windows (the engine of anti-evaporation) integrate
   to zero against the modular weight: negativity is a zero-sum reallocation.
3. `|I2|/I1` is a binary wedge-locality detector: 0 for boost-adapted data,
   1 for boost-blind data.

## Contents

| File | Purpose |
|---|---|
| `supplementary.tex` | Supplementary material: TikZ diagrams + numerical figures |
| `appendix_b_verification.py` | Numerical checks for Appendix B (I2 vanishing, negativity windows, one-mode sum rule) |
| `generate_figures.py` | Generates Figures S1-S3 (matplotlib) into `figures/` |
| `Makefile` | Build system |

## Quickstart

```sh
make verify        # run the numerical checks (numpy)
make figures       # regenerate Figures S1-S3 (matplotlib)
make               # build supplementary.pdf (needs pdflatex + TikZ)
```

Requirements: Python >= 3.10 with `numpy` and `matplotlib`; a TeX
distribution with TikZ (TeX Live recommended).

## Verified numbers

From `make verify` (kappa = 1):

* `|I2|/I1 = 1.3e-15` for a boost-positive-frequency packet (prediction: 0),
  with the boost/affine cross-check agreeing to `2e-16`.
* `|I2|/I1 = 1.000000` for a real (boost-blind) packet: no wedge-local
  squeeze exists at any finite r.
* Squeezed flux at `(r, theta) = (0.5, pi)`: negative on 51% of the packet
  support, funded depth 21% of the positive part, total positive.
* One-mode sum rule `S_rel = beta * dE` verified to <= 2e-15 at three
  parameter points, matching `beta * w * sinh(r)^2 * coth(beta*w/2)`.

## Status and open problems

* (i) Squeeze-shape inequality: **closed** (exact vanishing theorem).
* (ii) dS = 0 for the quadratic generator: verified at the one-mode level;
  the type-III domain-theoretic lift remains open.
* (iii) Backreaction with sign-indefinite flux (teleological horizon):
  the open dynamical frontier, constrained by the zero-sum theorem.

## Papers
- https://doi.org/10.5281/zenodo.21303604

