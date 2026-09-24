"""
Varied-response ("Method 2") normalisation for interaction dials.

Background
----------
For a cross-section measurement the signal rate in each true bin is a free
template parameter, measured from data rather than taken from the generator.
The conventional way of building a systematic universe -- reweight the events
with the dial and re-histogram them in reco space -- therefore charges the
measurement for a quantity it does not use: the generator's prediction of the
signal *rate*. Only the dial's effect on efficiency, bin migration and
backgrounds is a genuine uncertainty on the extraction.

This module implements the alternative. For dial ``k`` at sigma shift ``z`` and
an event ``e`` whose true value falls in true bin ``j``, the weight becomes

    w~_e(z) = w_e(z) / rho_j(z),

    rho_j(z) = sum_{true signal in bin j} w(z) / sum_{true signal in bin j} 1 ,

where the sums run over the *full* truth sample (selected or not), which is why
the truth tree and its systematic tree are needed. A dial that changes only the
normalisation of true bin j gives w = rho_j, hence w~ = 1 and no uncertainty --
correct, because a free template parameter already absorbs that change. A dial
that changes efficiency or migration reweights the selected events differently
from the true-bin average, so w~ != 1 and the variation survives.

The division is applied to signal events only. Background rates are taken from
the simulation, not measured, so backgrounds keep their unmodified weights.

Notes
-----
Sigma throws are derived deterministically from the dial name, so the selected
sample and the truth sample see the *same* universe for a given dial without
having to thread random state between them. This also makes runs reproducible.
Because the draw differs from an un-seeded ``np.random.normal`` call, absolute
numbers may shift by the usual O(1/sqrt(2*nuniv)) Monte Carlo noise relative to
previously produced plots; the Method 1 / Method 2 comparison is unaffected
since both use identical throws.
"""

import hashlib
import re

import numpy as np

# Fixed sigma grid that the interpolation is carried out on.
_SIGMA_GRID_7 = np.linspace(-3.0, 3.0, 7)
# GENIE multisigma weights are stored in this order.
_SIGMA_ORDER_6 = np.array([-1.0, 1.0, -2.0, 2.0, -3.0, 3.0])

_IDENT = re.compile(r'[A-Za-z_][A-Za-z0-9_]*')


def universe_sigmas(dial_name, nuniv, seed=0) -> np.ndarray:
    """
    Deterministic sigma throws for a single dial.

    Parameters
    ----------
    dial_name : str
        Name of the systematic dial. Used to seed the generator so that the
        selected sample and the truth sample draw identical universes.
    nuniv : int
        Number of universes.
    seed : int, optional
        Global offset, allowing the whole ensemble to be re-thrown.

    Returns
    -------
    numpy.ndarray
        Array of shape (nuniv,) of standard-normal sigma values.
    """
    digest = hashlib.md5(f'{dial_name}:{seed}'.encode()).digest()
    return np.random.default_rng(
        int.from_bytes(digest[:8], 'little')).standard_normal(nuniv)


def interpolate(weights_array, sigmas) -> np.ndarray:
    """
    Evaluate per-event sigma-point weights at a set of sigma values.

    Parameters
    ----------
    weights_array : numpy.ndarray
        Array of shape (nevents, 6) or (nevents, 7) of weights at the stored
        sigma points.
    sigmas : numpy.ndarray
        Array of shape (nuniv,) of sigma values to evaluate at.

    Returns
    -------
    numpy.ndarray
        Array of shape (nevents, nuniv).
    """
    npoints = weights_array.shape[1]
    if npoints == 6:
        order = np.argsort(_SIGMA_ORDER_6)
        levels = _SIGMA_ORDER_6[order]
        w = np.asarray(weights_array, dtype=float)[:, order]
    elif npoints == 7:
        levels = _SIGMA_GRID_7
        w = np.asarray(weights_array, dtype=float)
    else:
        raise ValueError(
            f'Cannot interpolate weights with {npoints} sigma points; '
            'expected 6 or 7.')

    w = np.where(np.isfinite(w), w, 1.0)

    out = np.empty((w.shape[0], len(sigmas)), dtype=float)
    for u, z in enumerate(sigmas):
        z = float(np.clip(z, levels[0], levels[-1]))
        hi = int(np.searchsorted(levels, z, side='left'))
        hi = min(max(hi, 1), len(levels) - 1)
        lo = hi - 1
        span = levels[hi] - levels[lo]
        frac = 0.0 if span == 0 else (z - levels[lo]) / span
        out[:, u] = w[:, lo] * (1.0 - frac) + w[:, hi] * frac
    return out


def _mask_columns(expression, available) -> list:
    """Branch names referenced by a pandas ``eval`` expression."""
    available = set(available)
    return [tok for tok in set(_IDENT.findall(expression)) if tok in available]


class ResponseNormalizer:
    """
    Computes and applies the rho_j(z) normalisation described above.

    Attributes
    ----------
    _edges : numpy.ndarray
        Bin edges of the true binning. This must be the binning of the
        template parameters in the fit, since that is what defines which
        rate changes are absorbed by the extraction.
    _tidx : numpy.ndarray
        True bin index of every event in the truth sample.
    _sig : numpy.ndarray
        Boolean mask selecting true signal events in the truth sample.
    _T_nom : numpy.ndarray
        Nominal population of each true bin (the rho denominator).
    """

    def __init__(self, file_handle, truth_tree, truth_systematics_tree,
                 truth_variable, truth_bin_edges, signal_mask,
                 nuniv=1000, seed=0) -> None:
        """
        Parameters
        ----------
        file_handle : uproot.reading.ReadOnlyDirectory
            Directory holding both the truth tree and its systematic tree,
            i.e. the same handle the Sample uses (``events/<key>``).
        truth_tree : str
            Name of the tree holding the *full* truth sample.
        truth_systematics_tree : str
            Name of the systematic-weight tree for that truth sample. It must
            be row-aligned with ``truth_tree``.
        truth_variable : str
            Branch defining the true binning, e.g. ``true_dpT_lp_genie``.
        truth_bin_edges : list
            Edges of the true binning.
        signal_mask : str
            A pandas ``eval`` expression selecting true signal events.
        nuniv : int, optional
            Number of universes. Must match the value used by Systematic.
        seed : int, optional
            Global throw offset.
        """
        self._edges = np.asarray(truth_bin_edges, dtype=float)
        self._nbins = len(self._edges) - 1
        self._nuniv = nuniv
        self._seed = seed
        self._truth_variable = truth_variable
        self._signal_mask = signal_mask
        self._cache = dict()

        tree = file_handle[truth_tree]
        self._syst = file_handle[truth_systematics_tree]

        # A mismatch here means the systematic tree was filled from a
        # different sample than the truth tree -- silently fatal, so refuse.
        if tree.num_entries != self._syst.num_entries:
            raise ValueError(
                f"Truth tree '{truth_tree}' has {tree.num_entries} entries but "
                f"'{truth_systematics_tree}' has {self._syst.num_entries}. "
                'These trees must be row-aligned for the varied-response '
                'method; check that the systematic tree was written from the '
                'truth sample and not from the selected sample.')

        columns = sorted(set([truth_variable]
                             + _mask_columns(signal_mask, tree.keys())))
        df = tree.arrays(columns, library='pd')

        self._sig = df.eval(signal_mask).to_numpy(dtype=bool)
        values = df[truth_variable].to_numpy()
        self._tidx = np.clip(np.digitize(values, self._edges) - 1,
                             0, self._nbins - 1)
        self._tidx[~np.isfinite(values)] = 0

        self._T_nom = np.bincount(self._tidx[self._sig],
                                  minlength=self._nbins).astype(float)

        empty = np.flatnonzero(self._T_nom == 0)
        if empty.size:
            print(f'[interaction_response] warning: true bins {list(empty)} '
                  'are empty in the truth sample; rho set to 1 there.')

    @property
    def seed(self):
        return self._seed

    @property
    def nuniv(self):
        return self._nuniv

    def rho(self, dial_name):
        """
        Per-true-bin population change for one dial.

        Returns
        -------
        numpy.ndarray or None
            Array of shape (ntruebins, nuniv), or None if the dial is not
            present in the truth systematic tree.
        """
        if dial_name in self._cache:
            return self._cache[dial_name]

        if dial_name not in self._syst:
            print(f'[interaction_response] warning: dial `{dial_name}` absent '
                  'from the truth systematic tree; left on the truth method.')
            self._cache[dial_name] = None
            return None

        raw = np.stack(self._syst[dial_name].array(library='np'))
        if raw.shape[1] not in (6, 7):
            self._cache[dial_name] = None
            return None

        sigmas = universe_sigmas(dial_name, self._nuniv, self._seed)
        weights = interpolate(raw, sigmas)

        varied = np.zeros((self._nbins, self._nuniv), dtype=float)
        np.add.at(varied, self._tidx[self._sig], weights[self._sig, :])

        denom = self._T_nom[:, np.newaxis]
        rho = np.divide(varied, denom, out=np.ones_like(varied), where=denom > 0)
        # A non-positive rho would flip the sign of the weight; leave those
        # events unscaled rather than producing nonsense.
        rho[~np.isfinite(rho) | (rho <= 0)] = 1.0

        self._cache[dial_name] = rho
        return rho

    def apply(self, universe_weights, sample, dial_name) -> np.ndarray:
        """
        Divide the signal events' universe weights by rho.

        Parameters
        ----------
        universe_weights : numpy.ndarray
            Array of shape (nevents, nuniv) for the selected sample.
        sample : Sample
            The parent Sample, used for the true binning and signal mask.
        dial_name : str
            Name of the dial being processed.

        Returns
        -------
        numpy.ndarray
            The rescaled universe weights.
        """
        rho = self.rho(dial_name)
        if rho is None:
            return universe_weights

        data = sample._data
        if self._truth_variable not in data.columns:
            raise ValueError(
                f'Truth binning variable `{self._truth_variable}` not present '
                f'in sample `{sample._name}`.')

        if rho.shape[1] != universe_weights.shape[1]:
            raise ValueError(
                f'rho was built with {rho.shape[1]} universes but the sample '
                f'was processed with {universe_weights.shape[1]}. Set '
                '`nuniv` in the interaction_method block to the same value '
                'used by Systematic.process.')

        signal = data.eval(self._signal_mask).to_numpy(dtype=bool)
        if signal.shape[0] != universe_weights.shape[0]:
            raise ValueError(
                f'Sample `{sample._name}` has {signal.shape[0]} rows but the '
                f'universe weights have {universe_weights.shape[0]}. The '
                'systematic weights and the sample data are out of step.')

        values = data[self._truth_variable].to_numpy()
        tidx = np.clip(np.digitize(values, self._edges) - 1, 0, self._nbins - 1)
        tidx[~np.isfinite(values)] = 0

        scale = np.ones_like(universe_weights)
        scale[signal, :] = rho[tidx[signal], :]
        return universe_weights / scale
