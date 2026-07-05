# %% [markdown]
# # Cathode Bending vs. proximity to cathode — reviewer comment #17 (Anne), take 3
#
# Confirmed from the previous exploration:
# - `events/full/selected` in the main analysis file is a real TTree (CV
#   sample) with `reco_vertex_x/y/z`, `reco_leading_muon_end_x/y/z`, etc.
# - The Cathode Bending variation (var09) in that same file is stored ONLY as
#   histograms (`var09_1D`, `var09_2D`) — no per-event position info.
# - The raw per-event Cathode Bending sample lives separately at:
#     /pnfs/icarus/persistent/users/spine/NuMI/Nu_det_sys/Cath_bend/*.root
#   per your sample config (`[[sample]] name = "var09"`).
#
# This script is exploratory-first again, because I don't know:
#   (a) whether that raw sample has the same branch names/tree structure as
#       `events/full/selected`,
#   (b) whether it needs xrootd instead of a direct POSIX path (pnfs/persistent
#       is usually POSIX-readable from FNAL interactive nodes, but not
#       necessarily from wherever you run this notebook),
#   (c) the real cathode-plane x-positions in your global coordinate
#       convention — I am not going to guess ICARUS geometry constants.
#
# Run Part 1 first, send me the output, then we fill in Part 2.

# %%
import glob
import uproot
import numpy as np
import matplotlib.pyplot as plt

MAIN_FILE = (
    "/Users/rvizarreta/Library/CloudStorage/GoogleDrive-rvizarreta14@gmail.com/"
    "My Drive/🏛 PhD Repository/🚀 Research/🤖 Experiments&Projects/ICARUS/"
    "ICARUS_CC0pi_Selection/data/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root"
)

VAR09_PATH_PATTERN = "/pnfs/icarus/persistent/users/spine/NuMI/Nu_det_sys/Cath_bend/*.root"

# %% [markdown]
# ## Part 1a — can we see the raw var09 files at all from here?
# If this list comes back empty, the path either isn't mounted in this
# environment (common — pnfs/persistent is usually only visible from FNAL
# interactive/gpvm nodes) or needs an xrootd URL instead of a plain path.
# Try the xrootd fallback in the commented block below if so.

# %%
var09_files = sorted(glob.glob(VAR09_PATH_PATTERN))
print(f"Found {len(var09_files)} files matching {VAR09_PATH_PATTERN}")
for fn in var09_files[:5]:
    print(" ", fn)

# If the glob above returns nothing, uncomment and adjust the xrootd form:
# import subprocess
# # e.g. list via `ifdh ls` or convert to root://fndcadoor.fnal.gov/ prefix:
# var09_files = [
#     "root://fndcadoor.fnal.gov/" + fn
#     for fn in glob.glob(VAR09_PATH_PATTERN.replace("/pnfs/", "/pnfs/"))
# ]

# %% [markdown]
# ## Part 1b — inspect one file's structure
# Don't assume it mirrors `events/full/selected` — check directly.

# %%
if var09_files:
    vf = uproot.open(var09_files[0])
    print("Top-level keys:")
    for k in vf.keys():
        print(" ", k)
else:
    print("No files found yet — resolve Part 1a before continuing.")

# %% [markdown]
# ## Part 1c — CORRECTION: this is a raw SBN CAF file, not the flat analysis
# ntuple format
#
# The Part 1b output shows `recTree`, `TotalPOT`, `TotalEvents`,
# `metadata/metatree` — that's the standard SBN "Common Analysis File" (CAF)
# structure, not the flattened `events/full/selected`-style tree from the
# main analysis file. `recTree` holds the full StandardRecord, which means
# branch names will be deeply nested CAF paths (e.g. things shaped like
# `rec.slc.reco.trk.start.x`), almost certainly NOT `reco_vertex_x`.
#
# `recTree` in a CAF can easily have several thousand branches, so don't
# print all of them — grep for likely candidates instead. Also, since this
# went through SPINE (not Pandora), the exact branch naming under `rec.slc.*`
# may be SPINE-specific; don't assume standard SBN CAF branch names either,
# just search and see what's actually there.

# %%
VAR09_TREE_NAME = "recTree"  # <-- corrected from Part 1b output

if var09_files:
    vtree = uproot.open(var09_files[0])[VAR09_TREE_NAME]
    var09_branches = vtree.keys()
    print(f"Total branches in {VAR09_TREE_NAME}: {len(var09_branches)}")

    def grep_caf_branches(patterns):
        hits = [b for b in var09_branches if any(p.lower() in b.lower() for p in patterns)]
        return sorted(hits)

    print("\nBranches matching position/vertex/track keywords:")
    for b in grep_caf_branches(["vtx", "vertex", "start", "end", "trk", "position", "cathode", "tpc", "muon"]):
        print(" ", b)

# %% [markdown]
# ## Part 1d — normalization inputs (TotalPOT / TotalEvents)
# These live at the top level of each file, not inside recTree, and are
# exactly what's needed to resolve the NORM_CV / NORM_VAR09 placeholder in
# Part 2c below — the Cathode Bending sample and the CV sample are separate
# productions and need to be POT-normalized against each other, not compared
# as raw counts.

# %%
if var09_files:
    total_pot = 0.0
    total_events = 0
    for fn in var09_files:
        ff = uproot.open(fn)
        total_pot += ff["TotalPOT"].member("fVal") if hasattr(ff["TotalPOT"], "member") else float(ff["TotalPOT"])
        total_events += int(ff["TotalEvents"].member("fVal") if hasattr(ff["TotalEvents"], "member") else float(ff["TotalEvents"]))
    print(f"Cathode Bending sample: TotalPOT summed over {len(var09_files)} files = {total_pot:.4e}")
    print(f"Cathode Bending sample: TotalEvents summed = {total_events}")
    print("\nNote: TotalPOT/TotalEvents may be stored as a TParameter or a simple")
    print("histogram/tree depending on CAF version — if the .member('fVal') access")
    print("above errors, print type(ff['TotalPOT']) and we'll adjust the read.")

# %% [markdown]
# # PART 2 — the actual proximity-to-cathode binning
#
# Fill these in from Part 1 output and your own geometry knowledge before
# running. I'm leaving placeholders rather than invented numbers.

# %%
# Branch names. The CV branches come from the flat MAIN_FILE ntuple
# (confirmed names). The var09 branches come from a raw CAF (recTree) and
# will almost certainly be DIFFERENT, nested names (e.g. something under
# `rec.slc...`) — do not assume these match CV. Fill in from the Part 1c
# grep output.
VTX_X_CV = "reco_vertex_x"
OBSERVABLE_CV = "reco_leading_muon_ke"

VTX_X_VAR09 = None       # <-- FILL IN from Part 1c grep output (CAF path, e.g. "rec.slc.???.x")
OBSERVABLE_VAR09 = None  # <-- FILL IN, the CAF equivalent of leading muon KE

assert VTX_X_VAR09 and OBSERVABLE_VAR09, (
    "Fill in the real CAF branch names from Part 1c before running Part 2 — "
    "recTree branch names will not match the flat ntuple's reco_* names."
)

# Cathode-plane x-position(s) in the same coordinate convention as the vertex
# branches above. ICARUS has one cathode per cryostat (2 cryostats) — this
# should be a list of those x-values, not a single number. Pull from your
# geometry service / existing containment-cut logic (reco_containment_cut
# and true_containment_cut already exist in the ntuple, so whatever defined
# those boundaries is the authoritative source — reuse it rather than
# re-deriving it here).
CATHODE_X_POSITIONS = []  # <-- FILL IN, e.g. [-210.29, 210.29] (placeholder shape only)

assert CATHODE_X_POSITIONS, "Fill in CATHODE_X_POSITIONS before running Part 2."

# %% [markdown]
# ## Part 2a — load CV and Cathode Bending samples

# %%
cv_tree = uproot.open(MAIN_FILE)["events/full/selected"]
cv_arrs = cv_tree.arrays([VTX_X_CV, OBSERVABLE_CV], library="np")
cv_x, cv_obs = cv_arrs[VTX_X_CV], cv_arrs[OBSERVABLE_CV]

# NOTE: with ~8000 files, uproot.concatenate over all of them at once may be
# slow / memory-heavy. Consider testing on var09_files[:50] first, then
# scaling up once the branch names and normalization are confirmed correct.
var09_arrs = uproot.concatenate(
    [f"{fn}:{VAR09_TREE_NAME}" for fn in var09_files],
    filter_name=[VTX_X_VAR09, OBSERVABLE_VAR09],
    library="np",
)
var_x, var_obs = var09_arrs[VTX_X_VAR09], var09_arrs[OBSERVABLE_VAR09]

print(f"CV sample: {len(cv_x)} events")
print(f"Cathode Bending sample: {len(var_x)} events")

# %% [markdown]
# ## Part 2b — distance to nearest cathode plane

# %%
def distance_to_cathode(x, cathode_positions):
    x = np.asarray(x)
    cathodes = np.asarray(cathode_positions)
    return np.min(np.abs(x[:, None] - cathodes[None, :]), axis=1)

cv_dist = distance_to_cathode(cv_x, CATHODE_X_POSITIONS)
var_dist = distance_to_cathode(var_x, CATHODE_X_POSITIONS)

print("Distance-to-cathode distribution (CV sample):")
counts, edges = np.histogram(cv_dist, bins=20)
for c, lo, hi in zip(counts, edges[:-1], edges[1:]):
    print(f"  [{lo:7.1f}, {hi:7.1f}) cm : {c}")

# %% [markdown]
# ## Part 2c — bin by distance, compute normalized Var/CV ratio
#
# IMPORTANT — normalization: the CV sample and the Cathode Bending sample are
# separate MC productions and are very unlikely to have the same POT/number
# of generated events. Part 1d gave you total_pot for the Cathode Bending
# sample from the CAF files directly — get the equivalent CV POT (there's a
# `POT` branch/object in the main file under `events/full/POT`, check its
# type the same way Part 1d did) and scale one sample to the other's POT
# before taking the ratio. Do NOT just take a raw-count ratio; that silently
# bakes in a normalization bias that would masquerade as a proximity effect.
NORM_CV = 1.0     # <-- e.g. 1.0 if CV is the reference
NORM_VAR09 = 1.0  # <-- e.g. cv_pot / var09_total_pot, once both are read

DIST_BINS = np.linspace(0, max(cv_dist.max(), var_dist.max()), 8)  # adjust after Part 2b

cv_yield, _ = np.histogram(cv_dist, bins=DIST_BINS)
var_yield, _ = np.histogram(var_dist, bins=DIST_BINS)

cv_yield_norm = cv_yield * NORM_CV
var_yield_norm = var_yield * NORM_VAR09

ratio = np.divide(
    var_yield_norm, cv_yield_norm,
    out=np.full_like(var_yield_norm, np.nan, dtype=float),
    where=cv_yield_norm > 0,
)
bin_centers = 0.5 * (DIST_BINS[:-1] + DIST_BINS[1:])

# %% [markdown]
# ## Part 2d — plot

# %%
fig, ax = plt.subplots(figsize=(6, 4))
ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
ax.plot(bin_centers, ratio, marker="o", color="#003087")
ax.set_xlabel("Distance from reconstructed vertex to nearest cathode plane [cm]")
ax.set_ylabel("Cathode Bending Var / CV (normalized)")
ax.set_title("Cathode Bending impact vs. proximity to cathode")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## Things to double check before trusting this
#
# 1. Vertex-only distance may not capture "crossing the cathode" for long
#    tracks that start far from it but cross it — consider also computing
#    distance using `reco_leading_muon_end_x` and taking whichever endpoint
#    is closer, or checking whether the track's x-range spans a cathode
#    plane at all (a boolean "crosses" flag), since Anne asked about both
#    proximity AND crossing.
# 2. Split contained vs. exiting (`reco_leading_muon_containment` /
#    `reco_containment_cut`) — the existing systematics tables already show
#    Cathode Bending differs between these two samples (7.8% vs 6.4%), so
#    mixing them here could wash out or exaggerate the proximity effect.
# 3. Confirm CV and var09 selections are the same signal definition/cuts —
#    if `reco_cut_type` or `true_nuisance_is_signal` selections differ
#    between the two productions, the comparison isn't apples-to-apples.