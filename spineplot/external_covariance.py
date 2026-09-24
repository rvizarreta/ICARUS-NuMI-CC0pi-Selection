"""
external_covariance.py

A Systematic whose covariance comes from a precomputed FRACTIONAL
covariance matrix instead of per-event universe weights read from the
sample's own ROOT file.

Motivation: the Geant4 hadron-reinteraction ("fate") universes only exist
in a separate MC production (icarus_..._ppfx_G4.root), whose events are
not the events of the sample being plotted, so they cannot be attached
event by event. The fractional covariance F[i][j] built from that
production for the GUNDAM fit (calc_geant4_covariance.py) is applied here
to the plotted sample's own bin contents y_i:

    Cov[i][j] = F[i][j] * y_i * y_j

which is exactly what the GUNDAM covariance ParameterSet does in the fit.

Conventions (matched to the rest of spineplot):
  * The covariance is built in the units of Systematic.process(), i.e. BEFORE
    the exposure scaling. Sample.set_weight() multiplies it by scale**2
    later, so y_i here is sum(ppfx_cv_weight) with no POT scale factor.
  * Bins are half-open [lo, hi), like BinParser and Systematic.process().
  * The GUNDAM matrices carry one extra overflow bin ([800, 99999]) that
    the plot variables do not have; it is trimmed (a sub-block of a
    covariance is the covariance of that sub-range). The bin edges are
    checked against the GUNDAM binning file when one is given, so a
    mismatch stops the run instead of silently mis-assigning bins.
  * A variable without an entry in `files` gets a zero covariance (no
    contribution), so recipes that combine systematics across all
    variables still work.
"""

import numpy as np
import uproot

from systematic import Systematic

# GUNDAM overflow bins are written as [800, 99999]
_OVERFLOW_MIN_UPPER_EDGE = 1.0e4


def read_tmatrixtsym(path):
    """Read the single '*_covariance' TMatrixTSym<double> in `path`."""
    with uproot.open(path) as f:
        keys = [k for k in f.keys(cycle=False, recursive=True)
                if k.endswith('_covariance')]
        if len(keys) != 1:
            raise ValueError(
                f"Expected exactly one '*_covariance' object in {path}, "
                f"found {keys}.")
        obj = f[keys[0]]
        n = int(obj.member('fNrows'))
        el = np.asarray(obj.member('fElements'), dtype=float)
    if el.size == n * (n + 1) // 2:
        # Packed upper triangle, row by row (verified against the
        # '*_correlation' object stored in the same files).
        m = np.zeros((n, n))
        k = 0
        for i in range(n):
            for j in range(i, n):
                m[i, j] = m[j, i] = el[k]
                k += 1
    elif el.size == n * n:
        m = el.reshape(n, n)
    else:
        raise ValueError(
            f"{path}: {el.size} elements is neither {n*(n+1)//2} (packed "
            f"symmetric) nor {n*n} (full) for a {n}x{n} matrix.")
    if not np.allclose(m, m.T):
        raise ValueError(f"{path}: covariance matrix is not symmetric.")
    return m


def read_bin_edges(path):
    """Edges of a GUNDAM 1D float binning file ('variables: x x' + rows)."""
    lo, hi = [], []
    with open(path) as f:
        lines = [l for l in f.read().splitlines() if l.strip()]
    if not lines[0].startswith('variables:'):
        raise ValueError(f"{path}: first line must start with 'variables:'.")
    for line in lines[1:]:
        a, b = line.split()
        lo.append(float(a))
        hi.append(float(b))
    if not np.allclose(lo[1:], hi[:-1]):
        raise ValueError(f"{path}: bins are not contiguous.")
    return np.array([lo[0]] + hi)


class ExternalFractionalCovariance(Systematic):
    """
    Parameters
    ----------
    name : str
        Systematic name (e.g. 'G4_proton').
    files : dict
        {variable name: path to a ROOT file holding a fractional
        covariance TMatrixTSym}.
    binning_files : dict, optional
        {variable name: path to the GUNDAM binning txt the matrix was
        built with}. If given, its edges are checked against the
        variable's.
    label : str, optional
    weight_branch : str
        Per-event weight used to form y_i (default 'ppfx_cv_weight',
        the CV weight the fractional covariance was normalised with).
    """

    def __init__(self, name, files, binning_files=None, label=None,
                 weight_branch='ppfx_cv_weight'):
        super().__init__(name, None, label)
        self._files = dict(files)
        self._binning_files = dict(binning_files or {})
        self._weight_branch = weight_branch
        self._std = 0.0

    def process(self, sample, mask, nuniv=1000):
        self._covariances = dict()
        for var_name, variable in self._variables.items():
            if not variable._bin_edges:
                continue
            edges = np.asarray(list(variable._bin_edges.values())[0], float)
            nb = len(edges) - 1
            key = f'{self._name}_{var_name}'

            if var_name not in self._files:
                self._covariances[key] = np.zeros((nb, nb))
                continue

            frac = read_tmatrixtsym(self._files[var_name])

            # Bin-edge safety check, then trim the GUNDAM overflow bin.
            if var_name in self._binning_files:
                g_edges = read_bin_edges(self._binning_files[var_name])
                if len(g_edges) - 1 != frac.shape[0]:
                    raise ValueError(
                        f"[{self._name}/{var_name}] matrix is "
                        f"{frac.shape[0]}x{frac.shape[0]} but the binning "
                        f"file has {len(g_edges) - 1} bins.")
                if not np.allclose(g_edges[:nb + 1], edges, atol=1e-3):
                    raise ValueError(
                        f"[{self._name}/{var_name}] plot edges {edges} do not "
                        f"match GUNDAM edges {g_edges[:nb + 1]}.")
            if frac.shape[0] == nb + 1:
                if self._binning_files.get(var_name) is not None:
                    if g_edges[-1] < _OVERFLOW_MIN_UPPER_EDGE:
                        raise ValueError(
                            f"[{self._name}/{var_name}] extra GUNDAM bin is "
                            f"not an overflow bin (upper edge {g_edges[-1]}).")
                frac = frac[:nb, :nb]
            elif frac.shape[0] != nb:
                raise ValueError(
                    f"[{self._name}/{var_name}] matrix is "
                    f"{frac.shape[0]}x{frac.shape[0]}; expected {nb} (plot "
                    f"bins) or {nb + 1} (plus overflow).")

            # Bin contents of the sample being plotted (no POT scale;
            # Sample.set_weight applies scale**2 to the covariance).
            data = sample._data[variable._key].to_numpy(dtype=float)
            if self._weight_branch in sample._data.columns:
                w = sample._data[self._weight_branch].to_numpy(dtype=float)
            else:
                w = np.ones(len(data))
            if variable.mask is not None:
                vmask = sample._data.eval(variable.mask).to_numpy(dtype=bool)
                data, w = data[vmask], w[vmask]
            idx = np.digitize(data, edges) - 1
            ok = (idx >= 0) & (idx < nb)
            y = np.zeros(nb)
            np.add.at(y, idx[ok], w[ok])

            self._covariances[key] = frac * np.outer(y, y)

            tot = y.sum()
            self._std = (np.sqrt(y @ frac @ y) / tot) if tot > 0 else 0.0
            sig = np.sqrt(np.diag(frac))
            print(f"[{self._name}] {var_name}: {nb} bins, fractional sigma "
                  f"{100*sig.min():.2f}%..{100*sig.max():.2f}%, "
                  f"integrated {100*self._std:.2f}%")
