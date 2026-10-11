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

## Still open

* The orientation of the modular flow is adopted, not rederived from KMS on
  this horizon.
* Only the UU component of the s-wave dilaton equation is used; its other
  components are not checked.
* The seed amplitude in the validity window `kappa Delta u << ln S_N` is a
  parametric estimate.
* l >= 1, non-coherent states, and the quantum derivation for modified gravity.
