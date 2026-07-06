# PLAN — void-content audit: does "voids are not empty" require a rerun? (pre-registered)

*Drafted 2026-07-06 by the review session, in answer to the user's question: real voids contain
dark matter and diffuse gas (lensing: central δ_m ≈ −0.3…−0.5 for watershed voids, ≈ −0.8 in the
deepest large-void interiors; FRB dispersion / Lyα confirm the diffuse baryons) — the Milne
(empty, δ = −1) void of the base model does not exist in nature. This document audits every
committed verdict for empty-void dependence, concludes **no full-suite rerun is required**, and
specifies the two cheap additions that ARE warranted, plus consolidated housekeeping. This is
the next-wave work order; execution stops at the WP-B Stage-0 gate for user review.*

---

## 1. Dependency audit — which results assume empty voids?

| committed verdict | assumes empty voids? | under real void contents | action |
|---|---|---|---|
| P1 tracker fails (paper 1) | tests Wiltshire's *published* empty-void model as-is | stands by construction | none |
| R1 = RECONCILES, fitted χ² 1396.06 (paper 2) | no — kinematic identities; H_v is *derived*, not assumed | χ² numbers stand; but the derived contrast becomes a physicality question (§2) | new probe (§2) |
| R2-final = SHAPE-UNAVAILABLE (paper 2) | no — measures volume *fractions* via thresholds | stands | none |
| FORCED_FVOBS_FAILS_BIC_BAR (paper 2) | no — same kinematic machinery | stands | none |
| b_req = 8.417% (paper 3) | no — pure data-side ratio | stands | none |
| WP-H2′ b_pred = 2.4% FAILS (paper 3) | E_max uses the derived H_v ⇒ effectively the empty-void **maximum** | b_pred only *decreases* → FAILS strengthens | sensitivity row (§3) |
| WP-B three-phase plan | empty phases were the base | already amended (A1: CEILING/PHYSICAL) | run per amended plan |

**Conclusion (the user's question answered): no rerun.** Every verdict was computed at the
mechanism's mathematical best case (maximal void-wall contrast), which is precisely what a
refutation should test: if the best case fails, all physical cases fail a fortiori. Real void
contents move every number in the direction that *widens* the gaps. What is missing is the
quantification of that statement — two cheap computations:

## 2. New probe — the contrast budget (lives in paper 2: `src/probes/contrast_budget.py`)

**Question.** The fitted history *implies* a void-wall expansion contrast
(H_v − H_w)/H̄₀ ≈ 0.52 at z = 0 (committed: `modelV_probeR.json → derived_backreaction_V`,
with the dimensionless excess (H_v−H_w)/⟨H⟩ ≈ 0.37–0.59 across z). Note the empty-void
theoretical maximum in an EdS background is H_v/H̄ = 3/2, i.e. excess 0.5 — **the required
contrast sits essentially at the physical ceiling**. Can *real* voids, with lensing-measured
matter contents, supply it?

**Method (analytic/ODE, hours).** Evolve a spherical open-patch (top-hat underdensity) void in
the standard background for central matter deficits spanning the lensing band
δ_m ∈ {−0.3, −0.5, −0.8} (typical watershed centrals bracketed by the deepest interiors;
Clampitt & Jain 2015, DES stacked void lensing as the provenance anchors — cite, and carry a
±0.1 band). Output the available excess (H_v − H̄)/H̄ (z) per depth; compare against the required
curve with its Δχ²≤1 band.

**Pre-registered expectations and falsifier.** Expected: available ≈ 0.2–0.35 vs required
≈ 0.45–0.5 → a ×1.5–2.5 physicality shortfall *independent* of the availability (R2) and bias
(WP-H2′) failures — a third, independent strike, sealing R1's "reconciles" as unphysical.
Falsifier (report at full volume if hit): if the nonlinear solve shows real-void excess reaching
the required band within the lensing depths, the strike is withdrawn and the discussion says so.
Artifacts: `probes_out/contrast_budget.json` + adversarial `verify_contrast_budget.json`
(independent re-solve, different integrator). Paper 2 discussion gets the result; paper 3 §7
(WP-N) cites it.

## 3. Sensitivity row — physical-void b_pred (lives in paper 3)

Recompute E_max with H_v capped by §2's available band (instead of the derived/empty maximum):
E_max_physical = γ̄₀ H_v^phys / H_dress − 1, then b_pred_physical = E_max_physical · ⟨φ⟩_HF.
Expected ≈ 1.0–1.7% vs b_req = 8.417% — the FAILS verdict strengthens from 3.4σ toward ~4σ
(measurement-only). Append as a labeled row to `bpred_survey_averaged.json` (do not overwrite
the CEILING row — the ceiling/physical pairing is the point), and mirror one line in
`verify_bpred_survey.json`'s recompute. Same falsifier discipline as §2.

## 4. Housekeeping (consolidated from the review; do before new probes)

1. `bpred_survey_averaged.json`: fill the leftover template string
   `"verdict_robust": "FAILS|PARTIAL|RESOLVES"` → `"FAILS"` (measurement-only, directional).
2. Add the r_hom* impossibility line: solve E_max·⟨φ⟩_HF(r_hom*) = b_req for r_hom*
   (expected ≈ 200–250 h⁻¹Mpc); state that the 2M++ profile (|δ(<r)| < 0.02 by ≈ 137) excludes
   it, and restate the envelope check one-sidedly (entire band below b_req; the current
   `nsigma_env_hi = 0.31 → RESOLVES` cell double-counts the φ systematic at the envelope edge).
3. Replace the `σ_global = 0.01 ESTIMATE` placeholder with the computed H̄₀ scale error from the
   joint fit (one profile scan), in both b_req σ stacks.
4. Standardize the verify-artifact verdict key (`verdict_of_verification` vs `overall`) across
   all verify JSONs — one key, machine-readable.

## 5. Execution order and stop point

1. §4 housekeeping (one agent, minutes-to-hours).
2. §2 contrast budget + verify (paper 2).
3. §3 physical b_pred row + verify (paper 3).
4. WP-B Stage 0 per [`PLAN_WPB_threephase.md`](PLAN_WPB_threephase.md) (as amended by A1 —
   report Q_avail under both galaxy-bias and lensing-depth mappings).
5. **STOP at the Stage-0 go/no-go** and report — Stage 1 (the dynamical solve) starts only on
   user review of the gate.

Fleet discipline: ≤ 6 agents, phase-by-phase, adversarial verification before any number is
quoted, failures at full volume. No paper tex is touched in this wave.

## 6. What this wave cannot change

The Hubble-tension verdict (WP-H2′ FAILS) and the availability verdict (SHAPE-UNAVAILABLE) are
not in play — §§2–3 can only strengthen or (per falsifier) locally withdraw one supporting
strike. The only live question in this wave is WP-B Stage 0's Q budget, and its ceiling is the
dark-energy-free geometry claim, not the tension.
