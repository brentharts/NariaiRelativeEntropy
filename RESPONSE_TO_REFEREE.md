# Response to the PLB-D-26-03020 referee report (v9 -> v10)

`nariai_relative_entropy_v10.tex` replaces v9. The checkable claims of the
report were tested first (`python3 check_referee_points.py`; every check has a
control that is made to fail on purpose). All of them hold. The paper has been
reframed accordingly: it no longer claims to prove the linear-response
postulate.

| # | Referee point | Verified? | Change in v10 |
|---|---|---|---|
| 1 | The JT "theorem" assumes the linearized Einstein equations; it is the converse direction already in Dorau-Much; the 8 pi and the 1/4 are the same G; Sec. 4.4 is circular | Agreed (also reached independently from the first-law bookkeeping: alpha = kappa dA/dE = 8 pi is the first law) | New title and abstract. Sec. 4.4/4.5 rewritten as a *consistency check* of the converse direction. Claim that the postulate "becomes a theorem" removed. Ansatz marked as taken from DM Eq. (19). |
| 2 | Sign: retarded zero-mode-free solution gives dA = -4 S_rel (focusing); the + sign needs teleological data or a kernel admixture | **Yes** (`[2]`): retarded phi(0) = -(8 pi/Lambda) * modular integral < 0. Teleological: memory term (>0) plus modular term with the *same negative sign* | Eqs. (retarded), (memory) added; boundary-condition dependence stated; only the magnitude and modular weighting are claimed. |
| 3 | Timescale hierarchy is wrong; Nariai perturbations grow classically at the Hubble rate; the Lambda term in Eq. (46) has the wrong sign | **Yes** (`[1]`): sinh T, cosh T cos(theta), cosh T sin(theta) solve the equation only with -Lambda g phi; they grow as exp(sqrt(Lambda) tau) | Sec. 2.4 rewritten: classical growth at rate kappa, window kappa Delta u << ln S_N (estimate, flagged as such). Sign and 1/Lambda factor of Eq. (46) corrected. The "largest hierarchy in physics" claim is removed. |
| 4 | Two distinct bifurcate Killing horizons, not antipodal cuts; static face integrates over both; antipodal sum rule applies the identification to the excluded kernel mode and needs "signed" S_rel | Agreed | Horizon structure restated (Sec. 1, 2.2, 5.1). Caption: H_B^R is the *past* horizon of R. Static/null face caveat added. Antipodal balance removed. |
| 5 | Nariai is not the maximal-entropy horizon; de Sitter has A = 12 pi/Lambda; Nariai is not unique; de Sitter offers the same "arena" | **Yes** (`[3]`): cosmological area decreases 12 pi/Lambda -> 4 pi/Lambda, black-hole area increases 0 -> 4 pi/Lambda | Extremality and uniqueness claims removed; ledger equation removed; novelty relative to D'Angelo et al. stated modestly. |
| 6 | Branch-selection rule is the focusing theorem; positivity of S_rel adds nothing; applying DA = 4 S_rel to dynamical cuts is out of scope | Agreed | Sec. 5 reduced to what is true; "predictions" removed. |
| 7 | Overlap with DM; reliance on unrefereed companion; Sec. 6/8 interpretive | Agreed | Overlap stated explicitly (Sec. 3 intro, DM equation numbers). Collapse speculation and ledger interpretation removed. Companion note marked unrefereed and unused. |
| m | Eq. (36) is the modular quadratic form, not the one-particle norm; Eq. (25) invokes Delta S on type III; K_R is a quadratic form; orientation of modular flow; Sec. 3.2 item 3 self-contradiction; l-blind kernel; m > 0 for Euclidean vacuum; Clausius relation is generic; Eq. (49) uses S = A/4; epsilon notation; beta is not a light-crossing time; references | Agreed | Each corrected in place (see diff). Orientation is now stated as adopted, not derived. References: journal DOIs added, [20] corrected to Cadamuro-Froeb-Katsinis-Mandrysch, Nariai originals cited, Bousso 1998 added. |

## One addition beyond the report

Sec. 4.5 records, from the companion CICY code (`pyCICY/theories/nariai.py`),
that for constant-curvature f(R) and quadratic gravity the coupling is
`alpha/pi = 8/F`, `dS_rel/dA = F/4 = dS_Wald/dA`, with
`F_eff = 1 + 2 Lambda (alpha + 4 beta)` from the Wald tensor. So the 8 pi is
Einstein's, not universal.

## Follow-up on the four items left open in v10 (-> v11)

`python3 check_open_items.py` (about 30 s; every check has a control, and a
sign-flipped mutation of the l >= 1 check was confirmed to fail).

| Item left open | Result in v11 |
|---|---|
| Orientation of the modular flow was adopted, not derived | **Derived** from the Tomita-Takesaki strip condition for the scaling-limit kernel on U < 0: with the dilation U' -> e^{2 pi t} U' the continuation runs from just below to just above the pole ray (correct boundary values); the opposite orientation interchanges them. The commutator is non-zero (Im part also obtained from the delta' term, a second route), so it is not a convention. Rests on the positive-energy (-i0) Kay-Wald kernel. Not compared term by term with DM Eq. (4). |
| Only the UU component of the dilaton equation used | **All 2D components checked** (TT, T-theta, theta-theta) of the linearized 4D Einstein equation about dS2 x S2 equal the JT operator with -Lambda g phi and normalization -(2/r0), independent of the 2D metric perturbation; the +Lambda sign fails. The angular component is the separate equation for the metric perturbation. |
| Seed amplitude is an estimate | **Not derived**, but shown to enter only logarithmically: for eps0 ~ S_N^-p the window is p ln S_N (140 / 281 / 562 for p = 1/2, 1, 2 at the observed Lambda). A seed that is not a power of S_N is not controlled. |
| l >= 1 not covered | **Covered at linear order**: the Raychaudhuri equation holds generator by generator, and S_rel is local in the transverse variable, so dA = -4 S_rel for any harmonic content (checked for l = 0..3 against the factorized sum). The localization step itself is DM's. |
| Non-coherent states not covered | **Covered for every wedge-local unitary excitation** psi = U Omega_0, U in N_R: the relative modular operator satisfies Delta_{Omega_0,psi} = Delta_{Omega_0} exactly, so S_rel = <psi, K psi> with no entropy difference (coherent and squeezed states are examples). Not covered: states not of that form (entangled across the bifurcation surface, mixed), where S_rel differs from Delta K by an entropy difference. |

## Follow-up 2: can S_rel = dA/4 be justified independently of the Einstein equations? (-> v12)

`python3 check_normalization.py` (about 2 s; controls included; a mutation of
the Deser-Tekin mapping was confirmed to fail).

* **No-go.** dA = c S_rel has c of dimension length^2; S_rel = 2 pi dE/kappa contains no G.
  So c/l_P^2 = 4 *defines* Newton's constant. Horizon kinematics plus QFT on a fixed
  background cannot produce the number 4.
* **What is not a definition:** that the same constant governs the horizon entropy and
  the asymptotic mass. Tested in L = (R + a Ric^2 + b R^2)/16 pi G by two code-disjoint
  routes: Wald entropy at the horizon (F = 1 + 2 Lambda (a + 4b)) against the
  Deser-Tekin charge at infinity (PRD 67, 084009, Eq. 36; literature, not rederived).
  They agree on a 48-point grid, vanish at the same critical point, and are independent
  of the horizon (r_b or r_c) and of M. A Gauss-Bonnet term shifts S by the constant
  64 pi^2 and does not change dS/dA.
* **Result:** with G_N defined by E = M/G_N, S = A/4G_N, alpha = 8 pi G_N, dA = 4 G_N S_rel,
  and the first law holds exactly on both horizons. So alpha = 8 pi is not "Einstein's"
  but alpha = 8 pi G_N is universal in this family. (v11's wording is corrected.)
* **Not established:** the value of G_N; completeness of the family (cubic and higher,
  non-local terms, extra modes that change the Newtonian limit); independence from the
  linearized field equations of each theory.

## Still open

* The normalization S_rel = dA/4 is not justified independently of the Einstein equations.
* Localization beyond linear order; states outside the wedge-local class.
* The quantum side of the relation for modified gravity.
* The amplitude of the quantum seed of the growing mode.
* K_R = 2 pi int (-U) T_UU as a quadratic form on the states considered is assumed (Bisognano-Wichmann).
