import re

import numpy as np
import pandas as pd

class Variable:
    """
    A class designed to encapsulate the configuration of a single
    variable for the analysis.

    Attributes
    ----------
    _name : str
        The name of the variable.
    _key : str
        The key/name of the branch in the input TTree containing the
        variable data.
    _range : tuple
        The range of the variable.
    _nbins : int
        The number of bins for the variable.
    _xlabel : str
        The x-axis label for the variable.
    _mask : string
        A mask formula to apply to the variable.
    _bin_edges : numpy.ndarray
        The bin edges for the variable.
    _bin_centers : numpy.ndarray
        The bin centers for the variable.
    _bin_widths : numpy.ndarray
        The bin widths for the variable.
    _validity_check : dict
        A dictionary containing the validity check for the variable
        in each sample. The key is the sample name and the value is
        a boolean indicating whether the variable is present in the
        sample. This is intended to be checked before using the
        variable in the analysis.
    """
    def __init__(self, name, key, range, nbins,
                 binning_scheme='equal_width', xlabel=None,
                 xlabel_fontsize=12,
                 mask=None, custom_bins = None, xticks=None) -> None:
        """
        Initializes the Variable object with the given kwargs.

        Parameters
        ----------
        name : str
            The name of the variable.
        key : str
            The key/name of the branch in the input TTree containing the
            variable data.
        range : tuple
            The range of the variable.
        nbins : int
            The number of bins for the variable.
        binning_scheme : str
            The binning scheme for the variable. This can be either
            'equal_width', 'custom_bins', or 'equal_population'. The 
            default is 'equal_width,' which creates bins of equal width
            irrespective of the number of entries in each bin.
            Option 'custom_bins' requires entering bin edges tuple for 
            parameter custom_bins.
        xlabel : str
            The x-axis label for the variable.
        mask : string, optional
            A mask formula to apply to the variable. The default is None.
        custom_bins : tuple
            The customized bin edges (required if binning_scheme is
            'custom_bins').
        xticks : list, optional
            Explicit x-axis tick positions. Matplotlib's automatic tick
            locator chooses round intervals that need not land on the
            end of the axis range, so e.g. an angular variable on
            [0, 180] can end up with its last label at 150 even though
            the axis and the data do extend to 180. Setting this forces
            the ticks (and therefore the end-of-range label). The
            default is None, which keeps the automatic locator.

        Returns
        -------
        None
        """
        self._name = name
        self._key = key
        self._range = range
        self._nbins = nbins
        self._binning_scheme = binning_scheme
        self._xlabel = xlabel
        self._xlabel_fontsize = xlabel_fontsize
        self._mask = mask
        self._custom_bins = custom_bins
        self._xticks = xticks
        self._validity_check = {}
        self._bin_edges = {}
        self._bin_centers = {}
        self._bin_widths = {}

    def check_data(self, categories, sample_name, sample):
        """
        Provides some functionality to check the data for the variable.
        Specifically, this function checks that each sample does indeed
        have the key for the variable. Additionally, it is also useful
        to parse the categories and provide a separate binning for each
        category group in the case of equal-population binning.

        If the variable is not found in a sample, a flag is set for
        the variable instance to indicate that the variable is missing
        and should not be used. If the variable is used later in the 
        analysis, the missing variable flag will be checked and an
        exception will be raised.

        Parameters
        ----------
        categories : dict
            A dictionary containing the categories for the analysis.
            The key is the category enumeration and the value is the
            name of the group that the enumerated category belongs to.
        sample_name : str
            The name of the Sample to be checked for validity.
        sample : Sample
            The Sample to be checked for validity/availability of the
            variable.

        Returns
        -------
        None.
        """
        self._validity_check[sample_name] = (self._key in sample._data.keys())

        if all(self._validity_check.values()):
            groups = {v: [] for v in categories.values()}
            for k, v in sample.get_data([self._key], self._mask)[0].items():
                if k in categories.keys():
                    groups[categories[k]].append(v[0])

            for g, v in groups.items():
                if v and self._binning_scheme == 'equal_population':
                    all_entries = pd.concat(v)
                    range_mask = ((all_entries >= self._range[0]) & (all_entries <= self._range[1]))
                    self._bin_edges[g] = np.percentile(all_entries[range_mask], np.linspace(0, 100, self._nbins+1))
                    self._bin_centers[g] = 0.5*(self._bin_edges[g][1:] + self._bin_edges[g][:-1])
                    self._bin_widths[g] = self._bin_edges[g][1:] - self._bin_edges[g][:-1]
                elif self._custom_bins is not None and self._binning_scheme == 'custom_bins':
                    self._nbins = len(self._custom_bins) - 1
                    self._range = (self._custom_bins[0], self._custom_bins[-1])
                    self._bin_edges[g] = np.array(self._custom_bins)
                    self._bin_centers[g] = 0.5*(self._bin_edges[g][1:] + self._bin_edges[g][:-1])
                    self._bin_widths[g] = self._bin_edges[g][1:] - self._bin_edges[g][:-1]
                else:
                    self._bin_edges[g] = np.linspace(self._range[0], self._range[1], self._nbins+1)
                    self._bin_centers[g] = 0.5*(self._bin_edges[g][1:] + self._bin_edges[g][:-1])
                    self._bin_widths[g] = self._bin_edges[g][1:] - self._bin_edges[g][:-1]

    # Trailing "(Contained + Exiting $\mu^\pm$)"-type suffix of an x label.
    _CONTAINMENT_RE = re.compile(
        r'\s*\(\s*(Contained \+ Exiting|Contained|Exiting)\s*\$[^$]*\$\s*\)\s*$')
    # Leading "Reco $...$" / "True $...$" of an x label.
    _RECO_TRUE_RE = re.compile(r'^(Reco|True)\s+(\$[^$]+\$)(.*)$')

    @property
    def hist_xlabel(self):
        """
        The x label used on 1D data/MC histograms. It is `xlabel` with
        (i) the muon-containment suffix removed (it is drawn as a text
        box on the plot instead, see `containment_label`) and (ii) a
        leading "Reco $X$" / "True $X$" rewritten as "$X^{reco}$" /
        "$X^{true}$". Labels that do not follow these patterns are
        returned unchanged.
        """
        label = self._xlabel
        if not isinstance(label, str):
            return label
        label = self._CONTAINMENT_RE.sub('', label)
        m = self._RECO_TRUE_RE.match(label)
        if m:
            tag = m.group(1).lower()
            math = m.group(2)[:-1] + r'^{\mathrm{' + tag + r'}}$'
            label = math + m.group(3)
        return label

    @property
    def containment_label(self):
        """
        The muon-containment selection quoted in `xlabel`, as a
        mathtext string ready to be drawn on a plot, or None if the
        label does not quote one.
        """
        if not isinstance(self._xlabel, str):
            return None
        m = self._CONTAINMENT_RE.search(self._xlabel)
        if m is None:
            return None
        text = m.group(1).replace(' ', r'\ ')
        return r'$\mathbf{' + text + r'\ \mu^{\pm}}$'

    @property
    def containment_text(self):
        """
        The muon-containment selection quoted in `xlabel` as plain text
        with a mathtext mu^pm (e.g. 'Contained + Exiting $\\mu^{\\pm}$'),
        or None if the label does not quote one.
        """
        if not isinstance(self._xlabel, str):
            return None
        m = self._CONTAINMENT_RE.search(self._xlabel)
        if m is None:
            return None
        return m.group(1) + r' $\mu^{\pm}$'

    @property
    def mask(self):
        """
        Getter method for the mask attribute.

        Parameters
        ----------
        None.

        Returns
        -------
        string
            The mask formula for the variable.
        """
        return self._mask