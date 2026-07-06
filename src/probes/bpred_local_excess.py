#!/usr/bin/env python3
"""b_pred -- the DECISIVE computation (paper-3 WP-H, PLAN P2).

The local Hubble bias a wall observer SHOULD measure, from the free-history
model's OWN two-scale structure, with ZERO new parameters, generalizing paper-1's
expansion-variance / apparent-Hubble-variance window (Wiltshire 2009: 17-22% for
the tracker at f_v0 ~ 0.7-0.8).

WHAT PAPER 1 DID (the machinery being generalized).
  Paper 1 never had a *script* that computed 17-22% -- that window is the
  Wiltshire(2009) apparent-Hubble-variance prediction (local maxima ~75 km/s/Mpc
  vs dressed asymptote ~61.7; freshH0.py hard-codes WINDOW=(0.17,0.22) as an
  anchor). The *computable* machinery paper 1/2 DO have is the exact two-scale
  kinematic decomposition (NOTES K1-K3, modelv_probeR.derived_curves):
      K1  Hbar = f_w H_w + f_v H_v                     (volume-avg expansion)
      K2  fv'  = 3 f_v (1-f_v)(H_v - H_w)              (fraction evolution)
      dressed:  H = gamma_bar Hbar - dgamma_bar/dt     (Wiltshire09 Eq 27, R3)
  This script GENERALIZES that: it builds the apparent-Hubble-variance window from
  the two-scale structure directly, VALIDATES that the construction reproduces
  Wiltshire's 17-22% when fed the TRACKER, then applies it to the free history.

THE CONSTRUCTION (zero new parameters).
  A wall observer's apparent Hubble rate varies with scale between two limits:
    * asymptotic (above the ~100 h^-1Mpc homogeneity scale): the dressed
      volume-average rate  H_dress = gamma_bar Hbar - gamma_bar'  (= Hd(0) Hbar0,
      the solver's present dressed rate; the GLOBAL H0 the ladder should converge to);
    * void-scale maximum (at the ~30 h^-1Mpc dominant-void scale, below homogeneity):
      the beam samples the fast void expansion, dressed to the wall clock,
      H_void_app = gamma_bar H_v - gamma_bar'.
  The predicted MAXIMUM apparent-Hubble excess (the free-history analogue of
  Wiltshire's 17-22% window) is
      b_pred_max = (H_void_app - H_dress) / H_dress = gamma_bar (H_v - Hbar)/H_dress.
  All quantities are exact kinematic identities of the forced f_v(z) at z=0 -- no
  fitted parameter enters.

  Several alternative zero-parameter excess measures are ALSO reported for
  transparency (bare void-over-mean, full void-wall contrast, pure lapse boost),
  so the validated calibration against the tracker window is auditable.

Writes probes_out/bpred_local_excess.json.
"""
import os, sys, json
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.abspath(os.path.join(_HERE, ".."))
os.chdir(_SRC)
sys.path.insert(0, _HERE)
sys.path.insert(0, _SRC)
import modelv_theory as MV

OUT = os.path.abspath(os.path.join(_SRC, "..", "probes_out", "bpred_local_excess.json"))

Z_NODES = [0.0, 0.3, 0.7, 1.3, 2.33]
FV_NODES_V = [0.64013, 0.53112, 0.39578, 0.27945, 0.19359]     # Probe R free-history winner (LA)
FV_NODES_V0 = [0.38292, 0.25043, 0.12798, 0.07152, 0.01699]    # V0 no-lapse control

# Wiltshire (2009) apparent-Hubble-variance window paper 1 cites for the tracker.
WILT_WINDOW = (0.17, 0.22)
# The paper-3 target (convention-independent bare ratio, roadmap +8.4%).
B_REQ_ROADMAP = 0.0842


def two_scale_at_z0(sol, lapse):
    """Exact two-scale decomposition at z=0 (all rates in units of Hbar0).

    Reproduces modelv_probeR.derived_curves at the z=0 node, independently.
    Returns a dict with H_w, H_v, Hbar, gamma_bar, gamma_bar_dot, H_dress and the
    candidate excess measures.
    """
    # locate z=0 (grid ascending)
    z = sol.z
    tau = sol.tau
    fv = sol.fv
    # fv0 and present values
    fv0 = float(sol.fv0)
    tau0 = float(np.interp(0.0, z, tau))
    fv_at0 = float(np.interp(0.0, z, fv))
    # df_v/dtau at z=0: fv'(z) * dz/dtau, both on the grid
    dz_dtau = np.gradient(z, tau)
    # df_v/dz from the callable is not stored on sol; reconstruct from grid fv(z)
    dfv_dz = np.gradient(fv, z)
    fvp = float(np.interp(0.0, z, dfv_dz * dz_dtau))            # df_v/dtau (>0)
    one_m = max(1.0 - fv0, 1e-9)

    Hw = 2.0 / (3.0 * tau0)                                     # H_w/Hbar0
    dHvw = fvp / (3.0 * fv0 * one_m)                            # (H_v - H_w)/Hbar0
    Hbar = Hw + fvp / (3.0 * one_m)                            # <H>/Hbar0  (=(1-fv)Hw+fv Hv)
    Hv = Hw + dHvw                                             # H_v/Hbar0

    if lapse == "algebraic":
        gam = (2.0 + fv0) / 2.0
        gamp = fvp / 2.0                                        # dgamma_bar/dtau
    elif lapse == "none":
        gam = 1.0
        gamp = 0.0
    else:
        raise ValueError(lapse)

    Hdress = gam * Hbar - gamp                                 # = Hd(0), the dressed present rate
    Hvoid_app = gam * Hv - gamp                                # dressed void rate (void-scale max)

    # candidate zero-parameter excess measures (all at z=0)
    E_dress_void = (Hvoid_app - Hdress) / Hdress               # PRIMARY: apparent-Hubble window max
    E_bare_void = (Hv - Hbar) / Hbar                           # bare void-over-mean
    E_contrast = dHvw / Hbar                                   # (H_v-H_w)/<H> (void_expansion_excess)
    E_lapse = gam - 1.0                                        # pure lapse boost fv0/2

    Hd0_solver = float(np.interp(0.0, z, sol.Hd))
    return dict(
        fv0=fv0, tau0=tau0, fvp_dtau=fvp,
        Hw_over_Hbar0=Hw, Hv_over_Hbar0=Hv, Hbar_over_Hbar0=Hbar,
        Hv_minus_Hw_over_Hbar0=dHvw,
        gamma_bar=gam, gamma_bar_dot=gamp,
        Hdress_over_Hbar0=Hdress, Hd0_solver_check=Hd0_solver,
        Hvoid_app_over_Hbar0=Hvoid_app,
        E_dress_void_PRIMARY=E_dress_void,
        E_bare_void=E_bare_void,
        E_contrast_HvHw_over_Hbar=E_contrast,
        E_lapse_boost=E_lapse,
    )


def solve_free(nodes, lapse):
    fv = MV.fv_from_nodes(nodes, z_nodes=Z_NODES)
    return MV.modelv_solve(fv, lapse=lapse, Ngrid=30000)


def solve_tracker(fv0, lapse="algebraic"):
    return MV.modelv_solve(MV.tracker_fv_of_z(fv0), lapse=lapse, Ngrid=30000)


def main():
    out = {"probe": "bpred_local_excess",
           "purpose": "b_pred: predicted local Hubble bias (wall observer) from the free-history "
                      "two-scale structure, zero new parameters; validated against Wiltshire's "
                      "17-22% tracker window.",
           "reading": "KINEMATIC (forced f_v(z); integrability not enforced).",
           "wiltshire_window": list(WILT_WINDOW),
           "b_req_roadmap": B_REQ_ROADMAP}

    # ---- VALIDATION: reproduce the tracker apparent-Hubble-variance window ----
    trk = {}
    for fv0 in [0.695, 0.76, 0.80, 0.85]:
        sol = solve_tracker(fv0, "algebraic")
        d = two_scale_at_z0(sol, "algebraic")
        trk[f"fv0={fv0}"] = d
    out["tracker_validation"] = trk
    # is the PRIMARY measure in the 17-22% window for the tracker at f_v0~0.7-0.8?
    prim_trk = {k: v["E_dress_void_PRIMARY"] for k, v in trk.items()}
    in_win = {k: bool(WILT_WINDOW[0] <= p <= WILT_WINDOW[1]) for k, p in prim_trk.items()}
    out["tracker_validation_primary"] = prim_trk
    out["tracker_validation_in_window"] = in_win

    # ---- b_pred for the FREE HISTORY (LA primary; V0 no-lapse control) ----
    solV = solve_free(FV_NODES_V, "algebraic")
    dV = two_scale_at_z0(solV, "algebraic")
    out["free_history_LA"] = dV

    solV0 = solve_free(FV_NODES_V0, "none")
    dV0 = two_scale_at_z0(solV0, "none")
    out["free_history_V0_nolapse"] = dV0

    # ---- validation-gate verdict: does the two-scale machinery reproduce 17-22%? ----
    # Wiltshire's window is the MAXIMUM apparent H0 ~72-75 vs dressed ~61.7. If the smooth
    # two-scale dressed-void-excess reproduced it, the tracker at f_v0~0.7-0.8 would land in
    # [0.17,0.22]. It does NOT (it gives ~0.09-0.12 and DECREASES with f_v0) -> the window is
    # a Wiltshire(2009) SPATIAL-variance literature result, not this codebase's machinery.
    tracker_gate_pass = any(in_win.values())
    out["window_machinery_status"] = dict(
        tracker_validation_gate_pass=bool(tracker_gate_pass),
        finding=("The paper-1 17-22% apparent-Hubble-variance window is a Wiltshire(2009) "
                 "LITERATURE CITATION (hard-coded WINDOW=(0.17,0.22) in freshH0.py); it is "
                 "never computed in the paper-1/2 codebase. The one computable expansion-"
                 "variance object present -- the exact two-scale kinematic decomposition "
                 "(NOTES K1-K3 / modelv_probeR.derived_curves) -- does NOT reproduce it: fed "
                 "the tracker it yields ~9-12%% at f_v0~0.7-0.8 (decreasing with f_v0), never "
                 "in [0.17,0.22]. The smooth volume-average carries no void-scale apparent-H0 "
                 "bump; that bump is a SPATIAL-inhomogeneity effect (radial apparent-H0 profile "
                 "/ observer environment) that both papers declare FUTURE / OUT-OF-SCOPE "
                 "(paper-1 tex decisive-test iii; PLAN P4 flow-correction re-derivation)."),
    )

    # ---- the two-scale BOUND (not a validated survey-averaged prediction) ----
    b_pred_max = dV["E_dress_void_PRIMARY"]
    # sigma per PLAN P3: anchored-scale (+) global-fit (+) lapse-reading spread, in quadrature.
    sig_anch = 0.7989638823689376 / 58.32194048848205          # phaseF free_fixed scale_err_sym/Hbar0
    sig_glob = 0.010                                            # global-fit Hbar0 rel err (SN+BAO+CMB, est.)
    # lapse-reading spread: LA (algebraic) vs V0 (no-lapse); LB (rate-ratio) NOT computed.
    lapse_hi = dV0["E_dress_void_PRIMARY"]                      # V0 no-lapse
    lapse_lo = b_pred_max                                       # LA
    sig_lapse = abs(lapse_hi - lapse_lo) / 2.0
    sigma = float(np.sqrt(sig_anch**2 + sig_glob**2 + sig_lapse**2))
    diff = float(b_pred_max - B_REQ_ROADMAP)
    nsig = float(abs(diff) / sigma)

    out["b_pred"] = dict(
        b_pred_max_apparent_bound_LA=b_pred_max,
        interpretation=("MAXIMUM apparent-Hubble excess at the dominant-void scale (z->0) from "
                        "the smooth two-scale structure -- a BOUND, not the survey-averaged "
                        "prediction the P3 equality test needs. The SH0ES ladder samples mostly "
                        "ABOVE the void scale, so the survey-averaged bias is a dilution of this "
                        "maximum toward b_req; computing that dilution needs the radial "
                        "apparent-H0 decay profile (missing). This BOUND is NOT validated: it "
                        "fails the tracker 17-22%% gate above."),
        b_pred_bare_void_LA=dV["E_bare_void"],
        b_pred_contrast_LA=dV["E_contrast_HvHw_over_Hbar"],
        b_pred_lapse_boost_LA=dV["E_lapse_boost"],
        b_pred_max_V0_nolapse=dV0["E_dress_void_PRIMARY"],
        b_req_roadmap=B_REQ_ROADMAP,
        window_bound_contains_b_req=bool(0.0 <= B_REQ_ROADMAP <= b_pred_max),
        sigma_components=dict(anchored_scale=float(sig_anch), global_fit=float(sig_glob),
                              lapse_spread_LA_V0=float(sig_lapse)),
        sigma_quadrature=sigma,
        diff_bound_minus_req=diff,
        nsigma_bound_vs_req=nsig,
        note_nsigma=("The LA bound sits %.2f sigma from b_req, but the 'agreement' is an "
                     "artifact of the enormous unresolved lapse spread (LA=%.3f vs V0=%.3f); "
                     "and the bound is the void-scale MAXIMUM, not the survey-averaged value. "
                     "Not a RESOLVES." % (nsig, lapse_lo, lapse_hi)),
    )

    out["verdict"] = dict(
        verdict="BLOCKED_needs_machinery",
        reason=("The decisive survey-averaged b_pred (the point prediction for the P3 equality "
                "test) is not computable with existing machinery. (1) Paper-1's 17-22%% window "
                "is a literature citation, never a codebase computation. (2) The available "
                "smooth two-scale expansion-variance machinery fails the tracker validation "
                "gate (gives ~9-12%%, not 17-22%%), so its free-history value (LA 19.5%%, V0 "
                "43.5%%) is neither validated nor the survey-averaged quantity b_req=8.4%% must "
                "be compared to. (3) The missing piece -- the radial apparent-H0 / spatial "
                "expansion-variance profile that turns the void structure into a survey-weighted "
                "local bias -- is exactly what both papers declare future/out-of-scope."),
        what_is_established=("controls reproduce (LCDM 73.53+-1.02, tracker 73.06); b_req=8.42%% "
                             "reproduced exactly and is convention-independent; the two-scale "
                             "(H_v-H_w)/Hbar0=0.516 at z=0 confirms the model has ample expansion "
                             "variance (a BOUND of 19.5%%-43.5%% CONTAINS b_req=8.4%%, i.e. the "
                             "mechanism is not shown incapable -- so NOT a FAILS)."),
    )

    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    # console summary
    print("=== TRACKER VALIDATION (primary E_dress_void vs Wiltshire 17-22%) ===")
    for k in trk:
        print(f"  {k}: b_pred_max={prim_trk[k]:.4f}  in_window={in_win[k]}  "
              f"(H_v/Hbar0={trk[k]['Hv_over_Hbar0']:.4f} Hdress={trk[k]['Hdress_over_Hbar0']:.4f})")
    print("=== FREE HISTORY (LA primary) ===")
    print(f"  fv0={dV['fv0']:.5f}  H_v/Hbar0={dV['Hv_over_Hbar0']:.5f}  "
          f"Hbar/Hbar0={dV['Hbar_over_Hbar0']:.5f}  Hdress/Hbar0={dV['Hdress_over_Hbar0']:.5f} "
          f"(solver Hd0={dV['Hd0_solver_check']:.5f})")
    print(f"  (H_v-H_w)/Hbar0={dV['Hv_minus_Hw_over_Hbar0']:.5f}  gamma_bar={dV['gamma_bar']:.5f} "
          f"gamma_bar_dot={dV['gamma_bar_dot']:.5f}")
    print(f"  b_pred_max (apparent window)   = {b_pred_max:.4f}")
    print(f"  b_pred bare_void               = {dV['E_bare_void']:.4f}")
    print(f"  b_pred contrast (Hv-Hw)/<H>    = {dV['E_contrast_HvHw_over_Hbar']:.4f}")
    print(f"  b_pred lapse boost             = {dV['E_lapse_boost']:.4f}")
    print(f"  b_pred V0 no-lapse max         = {dV0['E_dress_void_PRIMARY']:.4f}")
    print(f"  b_req (roadmap)                = {B_REQ_ROADMAP:.4f}")
    print(f"  tracker validation gate pass   = {out['window_machinery_status']['tracker_validation_gate_pass']}")
    print(f"  bound window contains b_req    = {out['b_pred']['window_bound_contains_b_req']}")
    print(f"  sigma (quadrature)             = {out['b_pred']['sigma_quadrature']:.4f}")
    print(f"  nsigma (bound vs b_req)        = {out['b_pred']['nsigma_bound_vs_req']:.2f}")
    print(f"  VERDICT                        = {out['verdict']['verdict']}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
