# Regenerate the flux and GENIE systematics tables (technote Tables 2 and 3)
# from the current systematics file, with explicitly defined columns:
#
#   delta_CV    = (mean_u T_u - N) / N : net fractional shift of the universe
#                 ensemble total yield relative to the unweighted CV. Nonzero
#                 only through the asymmetry of the dial response, since the
#                 weight at sigma = 0 is interpolated between the +-1 sigma
#                 knots rather than anchored at 1 (6-knot GENIE dials), or
#                 taken from the 0-sigma knot (7-knot flux dials).
#   delta_1sig  = std_u(T_u) / N : centered RMS spread of the universe totals.
#
# The non-centered RMS (what spineplot's print_sys reports) satisfies
# rms^2 = delta_CV^2 + delta_1sig^2 and is also printed for reference.
#
# Universe generation matches spineplot/medulla: 1000 sigma throws from
# N(0,1) per dial, weights interpolated at the sorted sigma knots. A fixed
# per-dial seed (crc32 of the dial name) makes the output reproducible.
#
# Group rows: delta_CV combines linearly (signed sum), delta_1sig in
# quadrature (dials assumed uncorrelated, as in the covariance treatment).

import os
import zlib

import numpy as np
import uproot

DATA_FILE = os.environ.get(
    "SPINE_DATA_FILE",
    "/Users/rvizarreta/Library/CloudStorage/GoogleDrive-rvizarreta14@gmail.com/"
    "My Drive/\U0001F3DB PhD Repository/\U0001F680 Research/\U0001F916 Experiments&Projects/"
    "ICARUS/ICARUS_CC0pi_Selection/data/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root",
)
TREE_DIR = "events/full"
NUNIV = 1000

# GENIE dial -> physics group (edit as needed; unmatched dials go to 'Other')
GENIE_GROUPS = {
    "CCQE": ["ZExpA1CCQE", "ZExpA2CCQE", "ZExpA3CCQE", "ZExpA4CCQE",
             "VecFFCCQEshape", "RPA_CCQE", "CoulombCCQE", "LowQ2Suppression"],
    "MEC": ["NormCCMEC", "NormNCMEC", "DecayAngMEC"],
    "NCEL": ["MaNCEL", "EtaNCEL"],
    "Resonant": ["MaCCRES", "MvCCRES", "MaNCRES", "MvNCRES", "NonRESBGvpCC1pi",
                 "NonRESBGvpCC2pi", "NonRESBGvpNC1pi", "NonRESBGvpNC2pi",
                 "NonRESBGvnCC1pi", "NonRESBGvnCC2pi", "NonRESBGvnNC1pi",
                 "NonRESBGvnNC2pi", "NonRESBGvbarpCC1pi", "NonRESBGvbarpCC2pi",
                 "NonRESBGvbarpNC1pi", "NonRESBGvbarpNC2pi",
                 "NonRESBGvbarnCC1pi", "NonRESBGvbarnCC2pi",
                 "NonRESBGvbarnNC1pi", "NonRESBGvbarnNC2pi",
                 "RDecBR1gamma", "RDecBR1eta", "Theta_Delta2Npi",
                 "ThetaDelta2NRad", "CC1piTPi"],
    "DIS/COH": ["NormCCCOH", "NormNCCOH", "AhtBY", "BhtBY", "CV1uBY", "CV2uBY"],
    "FSI": ["MFP_pi", "MFP_N", "FrCEx_pi", "FrInel_pi", "FrAbs_pi", "FrPiProd_pi",
            "FrCEx_N", "FrInel_N", "FrAbs_N", "FrPiProd_N"],
}


def dial_stats(tree, branch, nuniv=NUNIV):
    """Universe totals for one dial; returns (delta_cv, delta_1sig, rms) in %."""
    w = np.stack(tree[branch].array(library="np")).astype(float)
    n_events, n_knots = w.shape
    if n_knots == 1:
        # Single-alternative dial (e.g. DecayAngMEC): one weight per event
        # comparing the alternative model to the nominal. Treated as a
        # one-sided 1 sigma variation: delta_1sig = |yield shift|, no
        # central-value shift by construction.
        shift = 100 * (w[:, 0].sum() - n_events) / n_events
        return 0.0, abs(shift), abs(shift)
    if n_knots == 6:
        sigma_raw = np.array([-1, 1, -2, 2, -3, 3], dtype=float)
        order = np.argsort(sigma_raw)
        sigma_levels, w = sigma_raw[order], w[:, order]
    elif n_knots == 7:
        sigma_levels = np.linspace(-3, 3, 7)
    else:
        raise ValueError(f"{branch}: unexpected shape {w.shape}")

    rng = np.random.default_rng(zlib.crc32(branch.encode()))
    throws = rng.normal(0, 1, nuniv)
    # per-universe total yield: sum over events of the interpolated weight
    totals = np.einsum("eu->u", np.stack([np.interp(throws, sigma_levels, we)
                                          for we in w]))
    frac = (totals - n_events) / n_events
    d_cv = frac.mean()
    d_1s = frac.std()
    rms = np.sqrt((frac ** 2).mean())
    return 100 * d_cv, 100 * d_1s, 100 * rms


def fmt(x):
    return "$<$0.01" if abs(x) < 0.005 else f"{x:.2f}"


def run_family(tree, branches, groups, out_lines):
    results = {}
    for b in branches:
        short = b.split("multisigma_")[-1] if "multisigma" in b else \
                b.replace("hysyst_beam_", "").replace("hysyst_", "")
        results[short] = dial_stats(tree, b)
        print(f"  {short:>24}: d_CV = {results[short][0]:+6.2f}%  "
              f"d_1sig = {results[short][1]:5.2f}%  (rms {results[short][2]:5.2f}%)")

    used = set()
    for gname, dials in groups.items():
        present = [d for d in dials if d in results]
        if not present:
            continue
        used.update(present)
        g_cv = sum(results[d][0] for d in present)
        g_1s = np.sqrt(sum(results[d][1] ** 2 for d in present))
        out_lines.append(f"% --- {gname} ---")
        out_lines.append(rf"\rowcolor{{gray!20}}\multicolumn{{2}}{{l}}"
                         rf"{{\textbf{{\textit{{{gname}}}}}}} & "
                         rf"\textbf{{{fmt(g_cv)}}} & \textbf{{{fmt(g_1s)}}} \\")
        for d in present:
            out_lines.append(rf"\mbox{{{d}}} & & {fmt(results[d][0])} & "
                             rf"{fmt(results[d][1])} \\")
    leftover = [d for d in results if d not in used]
    if leftover:
        g_cv = sum(results[d][0] for d in leftover)
        g_1s = np.sqrt(sum(results[d][1] ** 2 for d in leftover))
        out_lines.append("% --- Other (unmatched dials, edit GENIE_GROUPS) ---")
        out_lines.append(rf"\rowcolor{{gray!20}}\multicolumn{{2}}{{l}}"
                         rf"{{\textbf{{\textit{{Other}}}}}} & "
                         rf"\textbf{{{fmt(g_cv)}}} & \textbf{{{fmt(g_1s)}}} \\")
        for d in leftover:
            out_lines.append(rf"\mbox{{{d}}} & & {fmt(results[d][0])} & "
                             rf"{fmt(results[d][1])} \\")
    return results


if __name__ == "__main__":
    f = uproot.open(DATA_FILE)
    out_dir = os.path.dirname(os.path.abspath(__file__))

    print("=== Flux dials (hysyst) ===")
    fx = f[f"{TREE_DIR}/selected_NuMIfluxsimTree"]
    fx_branches = sorted(k for k in fx.keys() if k.startswith("hysyst") and not k.endswith("_sigma"))
    flux_groups = {
        "Beam focusing": [b.replace("hysyst_beam_", "") for b in fx_branches
                          if b.startswith("hysyst_beam_")],
        "Hadron production (HPC)": [b.replace("hysyst_", "") for b in fx_branches
                                    if b.startswith("hysyst_hpc")],
    }
    lines_flux = []
    res_flux = run_family(fx, fx_branches, flux_groups, lines_flux)

    print("\n=== GENIE dials (multisigma) ===")
    ms = f[f"{TREE_DIR}/selected_multisigmaTree"]
    ms_branches = sorted(k for k in ms.keys()
                         if "multisigma_" in k and not k.endswith("_sigma"))
    lines_genie = []
    res_genie = run_family(ms, ms_branches, GENIE_GROUPS, lines_genie)

    with open(os.path.join(out_dir, "table_flux_systs.tex"), "w") as fh:
        fh.write("\n".join(lines_flux) + "\n")
    with open(os.path.join(out_dir, "table_genie_systs.tex"), "w") as fh:
        fh.write("\n".join(lines_genie) + "\n")
    print("\nWrote table_flux_systs.tex and table_genie_systs.tex")
