import numpy as np
from scipy.stats import binom
import pandas as pd
import matplotlib.pyplot as plt
import uproot

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary

class SpineEfficiency(SpineArtist):
    """
    A class designed to encapsulate the calculation of the efficiency
    of the selection as a function of the specified variable. The method
    employed is a Bayesian approach, where the true efficiency is assumed
    to be a binomial random variable with a uniform prior. The posterior
    distribution is then calculated using the binomial likelihood and the
    prior. 

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis. If None, the range is taken from the
        Variable object.
    _xtitle : str
        The title of the x-axis. If None, the title is taken from the
        Variable object.
    _variable : Variable
        The variable to calculate the efficiency with respect to.
    _samples : list
        A list of samples to use in the calculation of the efficiency.
    _categories : dict
        A dictionary mapping the category key to the category name.
    _cuts : dict
        A dictionary mapping the cut key to the cut label.
    _show_option : str
        The option to use when showing the artist.
    _npts : int
        The number of points to use when calculating the efficiency.
    _posteriors : dict
        A dictionary containing the posterior distributions for the
        efficiency calculation.
    _totals : dict
        A dictionary containing the total number of events in each bin
        of the variable.
    _successes : dict
        A dictionary containing the number of successful events in each
        bin of the variable.
    """

    # --- Hardcoded purity inputs --------------------------------------
    _PURITY_FILE = (
        "/Users/rvizarreta/Library/CloudStorage/GoogleDrive-"
        "rvizarreta14@gmail.com/My Drive/\U0001f3db PhD Repository/"
        "\U0001f680 Research/\U0001f916 Experiments&Projects/ICARUS/"
        "ICARUS_CC0pi_Selection/data/"
        "1muNp0pi_Nge1_uncontained_purity.root"
    )
    _PURITY_TREE_MAP = {
        "reco_fiducial_cut": "events/nominal/cut1_fiducial",
        "no_photons_cut": "events/nominal/cut3_no_photons",
        "no_electrons_cut": "events/nominal/cut4_no_electrons",
        "single_muon_cut": "events/nominal/cut5_single_muon",
        "no_pions_cut": "events/nominal/cut6_no_charged_pions",
        "multiple_protons_cut": "events/nominal/cut7_has_protons",
    }
    _PURITY_SIGNAL_CATEGORIES = (0, 1, 3, 4)
    def __init__(self, variable, categories, cuts, title,
                 xrange=None, xtitle=None, show_option='table',
                 npts=1e6):
        """
        Parameters
        ----------
        variable : Variable
            The variable to calculate the efficiency with respect to.
        categories : dict
            A dictionary mapping the category key to the category name.
        cuts : dict
            A dictionary mapping the cut key to the cut label.
        title : str
            The title of the artist.
        xrange : tuple, optional
            The range of the x-axis. If None, the range is taken from
            the Variable object. The default is None.
        xtitle : str, optional
            The title of the x-axis. If None, the title is taken from
            the Variable object. The default is None.
        show_option : str, optional
            The option to use when showing the artist. The default is
            'table.'
        npts : int, optional
            The number of points to use when calculating the efficiency.
            The default is 1e6.
        """
        super().__init__(title)
        self._xrange = xrange
        self._xtitle = xtitle
        self._variable = variable
        self._samples = list()
        self._categories = categories
        self._cuts = cuts
        self._show_option = show_option
        self._npts = int(npts)
        self._posteriors = dict()
        self._totals = dict()
        self._successes = dict()
        self._purity_results = None

    def draw(self, ax, show_option, percentage=True, show_seqeff=True,
             show_unseqeff=True, yrange=None, npts=1e6, style=None,
             logx=False, logy=False):
        """
        Draw the artist on the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        show_option : str, optional
            The option to use when showing the artist. The default is
            None. This is intended to be used in cases where the artist
            can be shown in different ways (e.g. 2D vs projection of 2D
            down to 1D).
        percentage : bool, optional
            A flag to indicate if the efficiency should be displayed
            as a percentage. The default is True.
        show_seqeff : bool, optional
            A flag to indicate if the sequential (cumulative)
            efficiency should be shown. The default is True.
        show_unseqeff : bool, optional
            A flag to indicate if the unsequential (single-cut)
            efficiency should be shown. The default is True.
        yrange : tuple, optional
            The range of the y-axis. The default is None.
        npts : int, optional
            The number of points to use when calculating the efficiency.
            The default is 1e6.
        style : Style, optional
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.
        logx : bool, optional
            A flag to indicate if the x-axis should be logarithmic.
            The default is False.
        logy : bool, optional
            A flag to indicate if the y-axis should be logarithmic.
            The default is False.

        Returns
        -------
        None.
        """
        ax.set_title(self._title)

        # Get a list of the groups, while respecting the order of the
        # groups as they are configured in the analysis block. It is
        # possible that duplicates exist (multiple categories within
        # the same group), so care must be taken to ensure that the
        # groups are unique AND in the correct order.
        groups = list()
        for category in self._categories.values():
            if category not in groups:
                groups.append(category)

        if show_option == 'table':
            # Lambda formatter to round the values to two decimal
            # places, display as percentage (if toggled), and add super
            # and subscripts for the error values.
            if percentage:
                formatter = lambda x,y,z: rf'${100*x:.2f}^{{\ +{100*y:.2f}}}_{{\ -{100*z:.2f}}}$'
                diff_key = 'Differential\nEfficiency [%]'
                cumu_key = 'Cumulative\nEfficiency [%]'
                purity_key = 'Purity [%]'
            else:
                formatter = lambda x,y,z: rf'${x:.2f}^{{\ +{y:.2e}}}_{{\ -{z:.2e}}}$'
                diff_key = 'Differential\nEfficiency'
                cumu_key = 'Cumulative\nEfficiency'
                purity_key = 'Purity'

            # Clear up the axis because we are going to draw a table
            # on it (no need for any other plot elements).
            ax.axis('off')

            purity_results = self._calculate_purity()

            def fmt_purity(cut_key):
                res = purity_results.get(cut_key)
                if res is None:
                    return r'$-$'
                cv, ms, ps = res
                return formatter(cv, ms, ps)

            # Create the table data.
            results = pd.DataFrame({r'   ': list(), r'Cut': list(), r'Efficiency [%]': list(), r'Cumulative [%]': list()})
            group_endpoint = dict()

            # Loop over the groups requested in the plot.
            for group in groups:
                # Extract the cv, msigma, and psigma values for the group.
                _, cv, msigma, psigma = self.reduce(group, significance=0.6827)

                seq = lambda x : [v for k, v in x.items() if 'unbinned_seq_' in k and 'unseq' not in k]
                unseq = lambda x : [v for k, v in x.items() if 'unbinned_unseq_' in k]

                entry = {   r'   ': [group,] + [r'' for _ in range(1, len(self._cuts))],
                            r'Cut': self._cuts.values(),
                            diff_key: [formatter(x,y,z) for x,y,z in zip(unseq(cv), unseq(msigma), unseq(psigma))],
                            cumu_key: [formatter(x,y,z) for x,y,z in zip(seq(cv), seq(msigma), seq(psigma))],
                            purity_key: [fmt_purity(cut_key) for cut_key in self._cuts.keys()]}
                results = pd.concat([results, pd.DataFrame(entry)])
                group_endpoint[group] = len(results)

            # Trim the table to only show the columns that are
            # requested by the user.
            cols = [r'   ', r'Cut',] if len(groups) > 1 else [r'Cut',]
            if show_unseqeff:
                cols.append(diff_key)
            if show_seqeff:
                cols.append(cumu_key)
            cols.append(purity_key)
            results = results[cols]

            # Rename "Cumulative" to "Efficiency" if it is the only
            # efficiency to show.
            if show_seqeff and not show_unseqeff:
                results.rename(columns={cumu_key : 'Efficiency [%]' if percentage else 'Efficiency'}, inplace=True)
            
            table_data = [results.columns.to_list()] + results.values.tolist()
            table = ax.table(cellText=table_data, colLabels=None, loc='center', cellLoc='left', edges='T')
            table.scale(1, 1.75)
            for i in range(2, len(table_data)):
                if i == len(table_data) - 1:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'B'
                else:
                    for j in range(len(table_data[i])):
                        table[i, j].visible_edges = 'open'
                    if i in group_endpoint.values():
                        table[i, 0].visible_edges = 'B'

            def calc_bbox_yext(obj):
                figure = plt.gcf()
                bbox = obj.get_window_extent(renderer=figure.canvas.get_renderer())
                p0, p1 = figure.transFigure.inverted().transform(bbox)
                return p1[1] - p0[1]

            scale = 2.75
            while calc_bbox_yext(table) > 0.92:
                table.scale(1, 1/scale)
                scale -= 0.05
                table.scale(1, scale)
            if scale < 2.0:
                print(f'Warning: Table with title `{self._title}` is too large to fit on the figure (scale = {scale:.2f}). Consider extending the figure vertically.')

            # Mark the POT and preliminary information on the plot.
            '''
            if style.mark_pot:
                mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=0.1)
            if style.mark_preliminary is not None:
                mark_preliminary(ax, style.mark_preliminary, vadj=0.1)
            '''

        elif show_option == 'differential':
            # Lambda formatter to round the values to two decimal
            # places, and display as percentage if toggled.
            if percentage:
                fmt = lambda x : 100*x
            else:
                fmt = lambda x : x

            # For clarity, this artist will only draw either the
            # sequential or the non-sequential efficiency. The
            # configuration file technically has a toggle for both,
            # but if both are requested, the artist will only draw the
            # sequential efficiency. The default is to draw the
            # sequential efficiency.
            if show_seqeff:
                key_base = 'binned_seq_'
            elif show_unseqeff:
                key_base = 'binned_unseq_'
            else:
                key_base = 'binned_seq_'

            # Note: if the user requests a differential plot, the
            # clarity of the plot is highly dependent on the number
            # of groups requested. If the number of groups is one,
            # it makes sense to denote the difference between cuts
            # with different colors. However, if the number of groups
            # is greater than one, it is better to color the groups
            # differently and denote the difference between cuts with
            # different marker styles. This can get a little messy,
            # but this gives the best "default" setting while leaving
            # the user the flexibility (and the responsibility) to 
            # ensure that the plot is clear and understandable.
            colors = ["green", "mediumblue", "brown", "goldenrod", "darkviolet", "dodgerblue"]
            if len(groups) == 1:
                for ci, (cut, cutname) in enumerate(self._cuts.items()):
                    _, cv, msigma, psigma = self.reduce(groups[0], significance=0.6827)
                    ax.errorbar(self._variable._bin_centers[groups[0]], fmt(cv[f'{key_base}{cut}']),
                                xerr=self._variable._bin_widths[groups[0]] / 2.0,
                                yerr=[  fmt(msigma[f'{key_base}{cut}']),
                                        fmt(psigma[f'{key_base}{cut}'])],
                                fmt='o', color=colors[ci], label=cutname, markersize=6, capsize=2, elinewidth=1, zorder=3)
            else:
                for gi, group in enumerate(groups):
                    for ci, (cut, cutname) in enumerate(self._cuts.items()):
                        _, cv, msigma, psigma = self.reduce(group, significance=0.6827)
                        ax.errorbar(self._variable._bin_centers[group], fmt(cv[f'{key_base}{cut}']),
                                    xerr=self._variable._bin_widths[group] / 2.0,
                                    yerr=[  fmt(msigma[f'{key_base}{cut}']),
                                            fmt(psigma[f'{key_base}{cut}'])],
                                    fmt=style.get_marker(ci), color=colors[gi],
                                    label=f'{group} : {cutname}', markersize=6, capsize=2, elinewidth=1, zorder=3)

            ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle, fontsize=12, weight='bold')
            ax.set_ylabel('Efficiency [%]' if percentage else 'Efficiency', fontsize=12, weight='bold')
            ax.set_xlim(self._variable._range if self._xrange is None else self._xrange)
            if yrange is not None:
                ax.set_ylim(yrange)

            # --- Background histogram of the true variable distribution ---
            # Uses the per-bin truth-level totals already accumulated for the
            # group during calculate(). Scaled so that the peak sits at ~50%
            # of the y-axis, with the bottom anchored at the current y-min.
            if len(groups) == 1 and groups[0] in self._totals:
                hist_group = groups[0]
                totals = np.asarray(self._totals[hist_group], dtype=float)
                if totals.max() > 0:
                    y_lo, y_hi = ax.get_ylim()
                    target_peak = y_lo + (y_hi - y_lo)
                    hist_scaled = y_lo + (totals / totals.max()) * (target_peak - y_lo)

                    bin_edges = self._variable._bin_edges[hist_group]
                    ax.bar(
                        bin_edges[:-1],
                        hist_scaled - y_lo,
                        width=np.diff(bin_edges),
                        bottom=y_lo,
                        align='edge',
                        color='grey',
                        alpha=0.15,
                        edgecolor='none',
                        zorder=0,
                    )
                    # Restore the y-limits in case bar() perturbed them.
                    ax.set_ylim(y_lo, y_hi)

            # Set the axis to be logarithmic if requested.
            if logx:
                # Modify the x-axis limits to ensure that the lower limit
                # is greater than zero. The lower edge needs to be at least
                # 3 orders of magnitude less than the maximum value in the
                # plot.
                xr = self._variable._range if self._xrange is None else self._xrange
                if xr == 0:    
                    xhigh_exporder = np.floor(np.log10(xr[1]))
                    xlow = xhigh_exporder - 3
                    ax.set_xlim(10**xlow, xr[1])
                ax.set_xscale('log')
            if logy:
                ax.set_yscale('log')

            ax.legend(fontsize=10)

            # Set tick mark size and tick label font size
            ax.tick_params(axis='both', which='major',
                           labelsize=12,  # Font size of tick labels
                           size=8,  # Length of major tick marks
                           width=2)  # Width/thickness of tick marks

            # Mark the POT and preliminary information on the plot.
            if style.scilimits:
                ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            if style.mark_pot:
                mark_pot(ax, self._exposure, style.mark_pot_horizontal)
            if style.mark_preliminary is not None:
                mark_preliminary(ax, style.mark_preliminary, hadj=0.035 if style.scilimits is not None else 0)

    def add_sample(self, sample, is_ordinate):
        """
        Add a sample to the efficiency calculation. The calculation of 
        the efficiency uses a Bayesian approach, where the true
        efficiency is assumed to be a binomial random variable with a
        uniform prior. The posterior distribution is then calculated
        using the binomial likelihood and the prior. The nature of this
        calculation means that it can be performed on a sample-by-sample
        (and category-by-category) basis, which is why the sample is
        processed individually.

        Parameters
        ----------
        sample : Sample
            The sample to add to the efficiency calculation.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        super().add_sample(sample, is_ordinate)
        self._samples.append(sample)

        # Calculate the efficiency for the sample
        self.calculate(sample)

    @staticmethod
    def multiply_posteriors(pos0, pos1):
        """
        Update the first posterior with the product of the two
        posterior distributions. This is a simple process, but it is
        important to normalize the posterior distribution after the
        multiplication. Additionally, the method is designed to handle
        both 1D and 2D posterior distributions (i.e. the method can
        handle both unbinned and binned efficiency calculations).

        Parameters
        ----------
        pos0 : np.ndarray
            The first posterior distribution.
        pos1 : np.ndarray
            The second posterior distribution.
        group0 : str
            The name of the first group, which is used to identify the
            posterior distribution if an exception is raised.
        group1 : str
            The name of the second group, which is used to identify the
            posterior distribution if an exception is raised.
        
        Returns
        -------
        pos : np.ndarray
            The product of the two posterior distributions.
        """
        pos = pos0*pos1
        if len(pos.shape) == 1:
            pos /= np.sum(pos)
        else:
            pos /= np.sum(pos, axis=-1)[:, np.newaxis]
        return pos

    def calculate(self, sample, significance=0.6827):
        """
        Calculate the efficiency of the selection on the sample as a
        function of the variable of interest. The method employed is a
        Bayesian approach, where the true efficiency is assumed to be a
        binomial random variable with a uniform prior. The posterior
        distribution is then calculated using the binomial likelihood
        and the prior. The nature of this calculation means that it can
        be performed on a sample-by-sample (and category-by-category)
        basis, which is why each category within the sample is processed
        individually.

        Parameters
        ----------
        sample : Sample
            The sample to calculate the efficiency for.
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        None.
        """
        # Create a linear space of efficiencies to calculate the
        # posterior distribution.
        efficiencies = np.linspace(0.0, 1, self._npts)

        # Get the data for the binning variable and the configured
        # cuts. The data is returned as a dictionary with the key
        # being the category and the value consisting of the variable
        # and the cuts.
        data, _ = sample.get_data([self._variable._key, *self._cuts.keys()], with_mask=self._variable.mask)
        for category, values in data.items():
            if category not in self._categories:
                continue
            if category not in self._posteriors:
                self._posteriors[category] = {f'seq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()}
                self._posteriors[category].update({f'unseq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()})
                

            # If the group does not already have an entry in the
            # posteriors dictionary, create one.
            if self._categories[category] not in self._posteriors:
                self._posteriors[self._categories[category]] = {f'unbinned_seq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()}
                self._posteriors[self._categories[category]].update({f'unbinned_unseq_{c}' : np.ones(efficiencies.shape) for c in self._cuts.keys()})

                self._totals[self._categories[category]] = np.zeros(self._variable._nbins)

                # If the show option is set to 'differential', create
                # the necessary dictionaries to store the binned
                # efficiencies.
                if self._show_option == 'differential':
                    
                    self._posteriors[self._categories[category]].update({f'binned_seq_{c}' : np.ones((self._variable._nbins, efficiencies.shape[0])) for c in self._cuts.keys()})
                    self._posteriors[self._categories[category]].update({f'binned_unseq_{c}' : np.ones((self._variable._nbins, efficiencies.shape[0])) for c in self._cuts.keys()})
                    
                    self._successes[self._categories[category]] = {f'binned_seq_{c}' : np.zeros(self._variable._nbins) for c in self._cuts.keys()}
                    self._successes[self._categories[category]].update({f'binned_unseq_{c}' : np.zeros(self._variable._nbins) for c in self._cuts.keys()})
            
            # The calculation of efficiency as a function of some the
            # variable of interest requires the binning of the variable
            # into a number of bins.
            bin_edges = self._variable._bin_edges[self._categories[category]]
            nbins = len(bin_edges) - 1

            # Calculate the success and total for each bin
            indices = np.digitize(values[0], bin_edges, right=False)
            self._totals[self._categories[category]] += np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
            for ci, (cut, cutname) in enumerate(self._cuts.items()):
                # Sequential cuts (unbinned)
                total = len(values[0])
                success = np.sum(np.all(values[1:ci+2], axis=0))
                self._posteriors[self._categories[category]][f'unbinned_seq_{cut}'] = SpineEfficiency.multiply_posteriors(self._posteriors[self._categories[category]][f'unbinned_seq_{cut}'], binom.pmf(success, total, efficiencies))
                
                # Non-sequential cuts (unbinned)
                success = np.sum(values[ci+1].to_numpy(bool))
                self._posteriors[self._categories[category]][f'unbinned_unseq_{cut}'] = SpineEfficiency.multiply_posteriors(self._posteriors[self._categories[category]][f'unbinned_unseq_{cut}'], binom.pmf(success, total, efficiencies))

                # If the show option is set to 'differential', calculate
                # the binned efficiencies.
                if self._show_option == 'differential':
                    # Sequential cuts (binned)
                    indices = np.digitize(values[0][np.all(values[1:ci+2], axis=0)], bin_edges, right=False)
                    success = np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
                    self._successes[self._categories[category]][f'binned_seq_{cut}'] += success
                    binomialpmf = [binom.pmf(self._successes[self._categories[category]][f'binned_seq_{cut}'][i],
                                             self._totals[self._categories[category]][i], efficiencies) for i in range(nbins)]
                    self._posteriors[self._categories[category]][f'binned_seq_{cut}'] = np.array(binomialpmf)

                    # Non-sequential cuts (binned)
                    indices = np.digitize(values[0][values[ci+1] == 1], bin_edges, right=False)
                    success = np.bincount(indices, minlength=len(bin_edges)+1)[1:-1]
                    self._successes[self._categories[category]][f'binned_unseq_{cut}'] += success
                    binomialpmf = [binom.pmf(self._successes[self._categories[category]][f'binned_unseq_{cut}'][i],
                                             self._totals[self._categories[category]][i], efficiencies) for i in range(nbins)]
                    self._posteriors[self._categories[category]][f'binned_unseq_{cut}'] = np.array(binomialpmf)

    def reduce(self, group, significance=0.6827):
        """
        Reduce the efficiency calculation to the final result (i.e.
        the central value and errors). This calculation sensibly
        treats the unbinned and binned efficiencies in the correct
        way: namely, the calculation preserves the shape of the
        unbinned and binned efficiencies. This is important for the
        usages downstream in the draw function.

        Parameters
        ----------
        group : str
            The key corresponding to the group of efficiencies to
            reduce.
        significance : float, optional
            The significance level to use when calculating the
            efficiency. The default is 0.6827.

        Returns
        -------
        final_posteriors : dict
            A dictionary containing the final efficiency posteriors
            for the requested group.
        cv : dict
            A dictionary containing the central value of the
            efficiency for the requested group.
        msigma : dict
            A dictionary containing the negative error on the
            efficiency for the requested group.
        psigma : dict
            A dictionary containing the positive error on the
            efficiency for the requested group.

        Raises
        ------
        ValueError
            If the group is not in the list of groups configured
            in the analysis block.
        """
        efficiencies = np.linspace(0.0, 1, self._npts)
        
        # If the group is not in the list of groups configured in
        # the analysis block, raise an exception. Otherwise,
        # extract the final posteriors for the group.
        if group not in self._posteriors.keys():
            raise ValueError(f"Group '{group}' not in the list of groups configured in the analysis block!")
        final_posteriors = self._posteriors[group]
        
        # Initialize the dictionaries to store the central value
        # and errors on the efficiency.
        cv = {}
        msigma = {}
        psigma = {}
        
        # The significance level is used to calculate the error on
        # the efficiency. The error is calculated by finding the
        # efficiency value that corresponds to the significance
        # level in the cumulative distribution of the posterior
        # distribution.
        sig = [0.5-significance/2.0, 0.5+significance/2.0]

        for key in final_posteriors.keys():
            if len(final_posteriors[key].shape) == 1:
                final_posteriors[key] /= np.sum(final_posteriors[key])
                cv[key] = efficiencies[int(np.argmax(final_posteriors[key]))]
                cumulative = np.cumsum(final_posteriors[key])
                msigma[key] = cv[key]-efficiencies[int(np.argmax(cumulative > sig[0]))]
                psigma[key] = efficiencies[int(np.argmax(cumulative > sig[1]))]-cv[key]

                # There is a case where posterior distribution
                # peaks at 0.0 or 1.0 due to 100% success or
                # failure. In this case, the error is set to 0.0,
                # which intuitively means that the error bar covers
                # only allowed values of the efficiency.
                if msigma[key] < 0:
                    msigma[key] = 0
                if psigma[key] < 0:
                    psigma[key] = 0

            else:
                final_posteriors[key] /= np.sum(final_posteriors[key], axis=1)[:, np.newaxis]
                cv[key] = efficiencies[np.argmax(final_posteriors[key], axis=1).astype(int)]
                cumulative = np.cumsum(final_posteriors[key], axis=1)
                msigma[key] = cv[key]-efficiencies[np.argmax(cumulative > sig[0], axis=1)]
                psigma[key] = efficiencies[np.argmax(cumulative > sig[1], axis=1)]-cv[key]

                # There is a case where posterior distribution
                # peaks at 0.0 or 1.0 due to 100% success or
                # failure. In this case, the error is set to 0.0,
                # which intuitively means that the error bar covers
                # only allowed values of the efficiency.    
                msigma[key][msigma[key] < 0] = 0
                psigma[key][psigma[key] < 0] = 0

        return final_posteriors, cv, msigma, psigma

    def _purity_posterior_stats(self, n_signal, n_total, significance=0.6827):
        """Bayesian posterior on purity, same grid + extraction as reduce()."""
        if n_total <= 0:
            return 0.0, 0.0, 0.0
        efficiencies = np.linspace(0.0, 1.0, self._npts)
        with np.errstate(divide='ignore', invalid='ignore'):
            log_post = (
                    n_signal * np.log(efficiencies)
                    + (n_total - n_signal) * np.log1p(-efficiencies)
            )
        log_post = np.where(np.isnan(log_post), -np.inf, log_post)
        log_post -= np.max(log_post)
        post = np.exp(log_post)
        post_sum = np.sum(post)
        if post_sum <= 0 or not np.isfinite(post_sum):
            return 0.0, 0.0, 0.0
        post /= post_sum

        cv = efficiencies[int(np.argmax(post))]
        cumulative = np.cumsum(post)
        sig_lo = 0.5 - significance / 2.0
        sig_hi = 0.5 + significance / 2.0
        msigma = cv - efficiencies[int(np.argmax(cumulative > sig_lo))]
        psigma = efficiencies[int(np.argmax(cumulative > sig_hi))] - cv
        if msigma < 0:
            msigma = 0.0
        if psigma < 0:
            psigma = 0.0
        return cv, msigma, psigma

    def _calculate_purity(self):
        """Read purity ROOT file and compute (cv, msigma, psigma) per cut."""
        if self._purity_results is not None:
            return self._purity_results

        results = {}
        try:
            f = uproot.open(self._PURITY_FILE)
        except (FileNotFoundError, OSError) as exc:
            print(f"Warning: purity file unreadable ({self._PURITY_FILE}): {exc}.")
            self._purity_results = {cut: None for cut in self._cuts.keys()}
            return self._purity_results

        try:
            for cut_key in self._cuts.keys():
                tree_path = self._PURITY_TREE_MAP.get(cut_key)
                if tree_path is None:
                    results[cut_key] = None
                    continue
                try:
                    cats = f[tree_path]["true_category"].array(library="np")
                except (KeyError, uproot.exceptions.KeyInFileError):
                    print(f"Warning: tree '{tree_path}' not found; skipping '{cut_key}'.")
                    results[cut_key] = None
                    continue

                n_total = int(len(cats))
                n_signal = int(np.isin(cats, self._PURITY_SIGNAL_CATEGORIES).sum())
                results[cut_key] = self._purity_posterior_stats(n_signal, n_total)
        finally:
            f.close()

        self._purity_results = results
        return self._purity_results