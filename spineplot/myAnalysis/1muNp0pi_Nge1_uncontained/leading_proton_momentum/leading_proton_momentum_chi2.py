# Data-MC chi2/ndf for the leading proton momentum -- reviewer question (London)
#
# London: "We should look at data-MC chi2/ndf comparisons for these plots
# (exiting and contained separately) to quantitatively check that the full
# model including all systematics can cover the distribution in data."
#
# Method: reuse the spineplot machinery so the covariance is IDENTICAL to the
# one behind the total-uncertainty band in the technote figures. The Analysis
# object is built from the same TOML; weights are set against the onbeam
# ordinate; each SpineSpectra1D artist accumulates the per-category stacked
# prediction and the data histogram; the 'MC_total' Systematic (flux + xsec +
# MC stat + detector, no G4) provides the covariance per variable. The data
# statistical uncertainty is added as a Poisson diagonal, diag(max(n_i, 1)).
#
#   chi2 = d^T C^-1 d,  d = data - prediction,
#   C = C_MC_total + diag(max(n_data, 1)),  ndf = number of bins.

import os
import sys

import numpy as np
from scipy import stats

workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))
from analysis import Analysis

data_file = os.environ.get(
    'SPINE_DATA_FILE',
        "/Users/rvizarreta/Library/CloudStorage/GoogleDrive-rvizarreta14@gmail.com/"
    "My Drive/\U0001F3DB PhD Repository/\U0001F680 Research/\U0001F916 Experiments&Projects/"
    "ICARUS/ICARUS_CC0pi_Selection/data/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root",
)

VARS = [
    'reco_leading_proton_momentum',
    'reco_leading_proton_momentum_contained',
    'reco_leading_proton_momentum_exiting',
]

ana = Analysis(os.path.join(workspace_root, 'leading_proton_momentum.toml'),
               data_file)

# Mimic Analysis.run() up to (but not including) figure creation: set the
# sample weights against the ordinate and let each artist accumulate its
# per-category histograms.
ordinate = ana._samples[ana._config['analysis']['ordinate_sample']]
for s in ana._samples.values():
    s.set_weight(target=ordinate)
for artist in ana._artists:
    for sample in ana._samples.values():
        artist.add_sample(sample, sample == ordinate)

# The 'MC_total' Systematic objects (one per sample that carries systematics)
# hold the covariance for every registered variable, already rescaled to the
# data exposure by set_weight.
systs = [s._systematics['MC_total'] for s in ana._samples.values()
         if 'MC_total' in getattr(s, '_systematics', {})]

print(f"\nSamples contributing an 'MC_total' covariance: "
      f"{[s._name for s in ana._samples.values() if 'MC_total' in getattr(s, '_systematics', {})]}")

done = set()
print(f"\n{'variable':>45} {'ndf':>4} {'chi2':>8} {'chi2/ndf':>9} {'p-value':>9}")
for artist in ana._artists:
    var = getattr(artist, '_variable', None)
    if var is None or var._name not in VARS or var._name in done:
        continue
    plotdata = getattr(artist, '_plotdata', None)
    if not plotdata or 'Data' not in plotdata:
        continue
    done.add(var._name)

    y_data = plotdata['Data']
    y_pred = np.sum([v for k, v in plotdata.items() if k != 'Data'], axis=0)

    cov = np.sum([s.get_covariance(var._name) for s in systs], axis=0)
    C = cov + np.diag(np.maximum(y_data, 1.0))

    d = y_data - y_pred
    chi2 = float(d @ np.linalg.solve(C, d))
    ndf = len(d)
    pval = stats.chi2.sf(chi2, ndf)
    print(f"{var._name:>45} {ndf:>4d} {chi2:>8.2f} {chi2 / ndf:>9.2f} {pval:>9.3f}")

    # Per-bin pulls for diagnostics: which bins drive the chi2
    sigma = np.sqrt(np.diag(C))
    pulls = d / sigma
    print(f"{'':>10} data  = {np.array2string(y_data, precision=1, max_line_width=200)}")
    print(f"{'':>10} pred  = {np.array2string(y_pred, precision=1, max_line_width=200)}")
    print(f"{'':>10} pulls = {np.array2string(pulls, precision=2, max_line_width=200)}")
