import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from spectra import SpineSpectra
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary, draw_error_boxes

class SpineSpectra2D(SpineSpectra):
    """
    A class designed to encapsulate a pair of variables' spectrum for
    an ensemble of samples. The class method add_sample() can be used
    to add a sample to the SpineSpectra2D.

    Attributes
    ----------
    _title : str
        The title of the artist. This will be placed at the top of the
        axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis for the spectrum. This is a tuple of
        the minimum and maximum values for the x-axis. If None, the
        range will be determined by Variable object assigned to the
        x-axis (show_option='2d') or set to (-1,1)
        (show_option='projection').
    _xtitle : str
        The label for the x-axis of the spectrum. If None, the label
        will be determined by the Variable object assigned to the
        x-axis (show_option='2d') or set to '(Y-X)/X'
        (show_option='projection').
    _yrange : tuple
        The range of the y-axis for the spectrum. This is a tuple of
        the minimum and maximum values for the y-axis. If None, the
        range will be determined by Variable object assigned to the
        y-axis (show_option='2d') or set to None
        (show_option='projection').
    _ytitle : str
        The label for the y-axis of the spectrum. If None, the label
        will be determined by the Variable object assigned to the
        y-axis (show_option='2d') or set to 'Entries'
        (show_option='projection').
    _variables : list
        The list of Variable objects for the spectrum.
    _categories : dict
        A dictionary of the categories for the spectrum. This serves as
        a map between the category label in the input TTree and the
        category label for the spectrum (and therefore what is shown
        in a single legend entry).
    _colors : dict
        A dictionary of the colors for the categories in the spectrum.
        This serves as a map between the category label for the
        spectrum (value in the `_categories` dictionary) and the color
        to use for the histogram. The color can be any valid matplotlib
        color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
    _plotdata : dict
        A dictionary of the data for the spectrum. This is a map between
        the category label for the spectrum and the histogram data for
        that category.
    """
    def __init__(self, variables, categories, colors, category_types,
                 title=None, xrange=None, xtitle=None, yrange=None,
                 ytitle=None, signal_categories=None) -> None:
        """
        Initializes the SpineSpectra2D object.

        Parameters
        ----------
        variables : list
            The list of Variable objects for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This serves
            as a map between the category label in the input TTree and
            the category label for the spectrum (and therefore what is
            shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the spectrum.
            This serves as a map between the category label for the
            spectrum (value in the `_categories` dictionary) and the color
            to use for the histogram. The color can be any valid matplotlib
            color string or a cycle indicator (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the spectrum.
            This serves as a map between the category label for the spectrum
            (value in the `_categories` dictionary) and the type of plot to
            use for the histogram. The type should be either 'histogram' or
            'scatter' to correspond to a stacked histogram or scatter plot,
            respectively.
        title : str, optional
            The title of the artist. This will be placed at the top of
            the axis assigned to the artist. The default is None.
        xrange : tuple, optional
            The range of the x-axis for the spectrum. This is a tuple of
            the minimum and maximum values for the x-axis. If None, the
            range will be determined by Variable object assigned to the
            x-axis (show_option='2d') or set to (-1,1)
            (show_option='projection'). The default is None.
        xtitle : str, optional
            The label for the x-axis of the spectrum. If None, the label
            will be determined by the Variable object assigned to the
            x-axis (show_option='2d') or set to '(Y-X)/X'
            (show_option='projection'). The default is None.
        yrange : tuple, optional
            The range of the y-axis for the spectrum. This is a tuple of
            the minimum and maximum values for the y-axis. If None, the
            range will be determined by Variable object assigned to the
            y-axis (show_option='2d') or set to None
            (show_option='projection'). The default is None.
        ytitle : str, optional
            The label for the y-axis of the spectrum. If None, the label
            will be determined by the Variable object assigned to the
            y-axis (show_option='2d') or set to 'Entries'
            (show_option='projection'). The default is None.
        signal_categories : list, optional
            The raw category values (as they appear in the category
            branch) that are to be treated as signal. Only used by
            show_option='smearing'. The response matrix in the forward
            model mu_j = sum_i R_ji s_i + b_j maps a true SIGNAL rate to
            reconstructed bins; backgrounds enter through b_j and have no
            signal truth bin, so including them dilutes the diagonal and
            the plotted matrix is then not the R the fit uses. If None,
            every category is summed, which reproduces the previous
            behaviour. The default is None.

        Returns
        -------
        None.
        """
        super().__init__(variables, categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
        self._category_types = category_types
        self._plotdata_diagonal = None
        self._binedges_diagonal = None
        self._signal_categories = signal_categories

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra2D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra2D object.
        is_ordinate : bool
            A flag to indicate if the sample is the ordinate sample.

        Returns
        -------
        None.
        """
        
        super().add_sample(sample, is_ordinate)

        if self._plotdata is None:
            self._plotdata = {}
            self._binedges = {}
        if self._plotdata_diagonal is None:
            self._plotdata_diagonal = {}
            self._binedges_diagonal = {}

        # Check if a mask is present for the variables. If so, we need
        # to combine the masks for the two variables.
        if self._variables[0].mask is not None and self._variables[1].mask is not None:
            joint_mask = f'{self._variables[0].mask} and {self._variables[1].mask}'
        elif self._variables[0].mask is not None:
            joint_mask = self._variables[0].mask
        elif self._variables[1].mask is not None:
            joint_mask = self._variables[1].mask
        else:
            joint_mask = None

        data, weights = sample.get_data([self._variables[0]._key, self._variables[1]._key], joint_mask)

        for category, values in data.items():
            if category not in self._categories.keys():
                continue
            xr = self._variables[0]._range if self._xrange is None else self._xrange
            yr = self._variables[1]._range if self._yrange is None else self._yrange
            x_bin_edges = list(self._variables[0]._bin_edges.values())[0] if self._variables[
                0]._bin_edges else np.linspace(xr[0], xr[1], self._variables[0]._nbins + 1)
            y_bin_edges = list(self._variables[1]._bin_edges.values())[0] if self._variables[
                1]._bin_edges else np.linspace(yr[0], yr[1], self._variables[1]._nbins + 1)
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros((len(x_bin_edges) - 1, len(y_bin_edges) - 1))
            h = np.histogram2d(values[0], values[1], bins=(x_bin_edges, y_bin_edges), weights=weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = (h[1], h[2])

            if self._categories[category] not in self._plotdata_diagonal:
                self._plotdata_diagonal[self._categories[category]] = np.zeros(self._variables[0]._nbins)
            diag = np.divide(values[1] - values[0], values[0])
            xr = (-1, 1) if self._xrange is None else self._xrange
            h = np.histogram(diag, bins=self._variables[0]._nbins, range=xr, weights=weights[category])
            self._plotdata_diagonal[self._categories[category]] += h[0]
            self._binedges_diagonal[self._categories[category]] = h[1]

    def draw(self, ax, style, show_option='2d', draw_identity=True,
             draw_colorbar=True, invert_stack_order=False,
             fit_type=None, logx=False, logy=False, logz=False,
             draw_stat_error=False, annotate_threshold=0.1,
             annotate_min_cell=0.022) -> None:
        """
        Plots the data for the SpineSpectra2D object.

        Parameters
        ----------
        ax : Axes
            The matplotlib Axes object to draw the plot on.
        style : Style
            The Style object to use for the plot.
        show_option : str
            The option to use for the plot. This can be one of a few
            options (default is '2d'):
                '2d'         - Draw a 2D histogram of the data.
                'projection' - Draw a projection of the data about the
                               diagonal.
        draw_identity : bool
            A flag to indicate if the identity line should be drawn on
            the plot. The default is True.
        draw_colorbar : bool
            A flag to indicate if a colorbar should be drawn on the plot.
            The default is True
        invert_stack_order : bool
            A flag to indicate if the stack order in the legend should
            be inverted. The default is False.
        fit_type : str
            The type of fit to perform on the data. The default is
            None, which will not perform any fit. The options are:
                'crystal_ball' - Perform a Crystal Ball fit on the data.
                'gaussian'     - Perform a Gaussian fit on the data.
        logx : bool
            A flag to indicate if the x-axis should be logarithmic.
            The default is False.
        logy : bool
            A flag to indicate if the y-axis should be logarithmic.
            The default is False.
        logz : bool
            A flag to indicate if the z-axis (colorbar) should be
            logarithmic. The default is False.
        draw_stat_error : bool
            A flag to indicate if the statistical error should be drawn
            on the plot. The default is False.
        annotate_threshold : float
            Only used by show_option='smearing'. Cells whose column
            normalized value (in percent) is at or below this are left
            unlabelled, so that the near-empty off-diagonal cells do not
            crowd the informative ones. The colour scale still shows
            them. Set to 0 to label every cell. The default is 0.1.
        annotate_min_cell : float
            Only used by show_option='smearing'. Cells whose width or
            height is below this fraction of the corresponding axis span
            are left unlabelled, because there is no room for the text.
            Needed for variable-width binnings where some bins are a
            small percentage of the axis. Set to 0 to disable. The
            default is 0.03.

        Returns
        -------
        None.
        """
        ax.set_title(self._title)
        
        if show_option == '2d' and self._plotdata is not None:
            values = np.sum([v for v in self._plotdata.values()], axis=0)
            binedges = self._binedges[list(self._plotdata.keys())[0]]
            
            # Find the minimum power of ten that is higher than the
            # maximum value in the plot. This will be used to set the
            # colorbar limits. The power of vmax is then given by the
            # maximum of this value and 3.
            max_power = np.floor(np.log10(np.max(values)))
            max_power = max([max_power, 2])
            ln = LogNorm(vmin=1, vmax=10**max_power)

            ax.imshow(values.T, extent=(binedges[0], binedges[-1], binedges[0], binedges[-1]),
                      aspect='auto', origin='lower', norm=ln if logz else None)
            ax.set_xlabel(self._variables[0]._xlabel if self._xtitle is None else self._xtitle, fontsize=12, weight='bold')
            ax.set_ylabel(self._variables[1]._xlabel, fontsize=12, weight='bold')
            ax.set_aspect('equal')

            # Set tick mark size and tick label font size
            ax.tick_params(axis='both', which='major',
                           labelsize=12,  # Font size of tick labels
                           size=8,  # Length of major tick marks
                           width=2)  # Width/thickness of tick marks
            
            # Draw the identity line. This must span the full range
            # of the plot, so we need to find the minimum and maximum
            # of the range for the plot.
            if draw_identity:
                min_range = min([binedges[0], binedges[-1]])
                max_range = max([binedges[0], binedges[-1]])
                ax.plot([min_range, max_range], [min_range, max_range], 'k--')

            # Draw the colorbar if requested. The color axis may also
            # be logarithmic if requested.
            if draw_colorbar:
                cbar = plt.colorbar(ax.images[0], ax=ax)
                cbar.set_label('Entries', fontsize=12, weight='bold')
                cbar.ax.tick_params(labelsize=12, size=6, width=1.5)

        if show_option == 'smearing' and self._plotdata is not None:
            # The smearing matrix is the response matrix R_ji of the forward model
            # mu_j = sum_i R_ji s_i + b_j, which is defined for SIGNAL only: a
            # background event has no signal truth bin to be unfolded to, and enters
            # the prediction through b_j instead. Summing every category dilutes the
            # diagonal and produces something that is not the R the fit uses.
            if self._signal_categories is not None:
                wanted = {lab for cat, lab in self._categories.items()
                          if any(float(cat) == float(s) for s in self._signal_categories)}
                if not wanted:
                    raise ValueError(
                        f'signal_categories={self._signal_categories} matched no entry '
                        f'in the category table {sorted(self._categories)}.')
                keys = [k for k in self._plotdata if k in wanted]
                if not keys:
                    raise ValueError(
                        f'signal_categories={self._signal_categories} resolved to labels '
                        f'{sorted(wanted)}, none of which are present in this figure.')
                values = np.sum([self._plotdata[k] for k in keys], axis=0)
                binedges = self._binedges[keys[0]]
            else:
                values = np.sum([v for v in self._plotdata.values()], axis=0)
                binedges = self._binedges[list(self._plotdata.keys())[0]]

            # Handle different binedges structures
            if isinstance(binedges, (list, tuple)) and len(binedges) == 2:
                x_edges = binedges[0]
                y_edges = binedges[1]
                extent = (x_edges[0], x_edges[-1], y_edges[0], y_edges[-1])
            else:
                x_edges = y_edges = binedges
                extent = (binedges[0], binedges[-1], binedges[0], binedges[-1])

            # Column normalize
            column_sums = np.sum(values, axis=1, keepdims=True)
            column_sums[column_sums == 0] = 1
            normalized_values = values / column_sums

            # ny corresponds to reco (y_edges), nx corresponds to truth (x_edges)
            ny, nx = normalized_values.T.shape
            x_bins = x_edges  # truth bins, nx+1 edges
            y_bins = y_edges  # reco bins, ny+1 edges

            # Create the plot using pcolormesh for correct variable-width bin handling
            im = ax.pcolormesh(x_edges, y_edges, normalized_values.T,
                               cmap='Reds', vmin=0, vmax=0.25)

            # Annotate each cell, subject to two suppressions. Both exist because a
            # label that cannot be read is worse than no label: it overlaps its
            # neighbours and obscures the colour underneath.
            #   1. VALUE: cells at or below annotate_threshold (percent) carry no
            #      information worth a number. With fine binning these near-empty
            #      off-diagonal cells vastly outnumber the informative ones.
            #      NOTE this can blank an ENTIRE row or column, not just scattered
            #      cells: the 2:1 median split routinely produces one reco bin much
            #      narrower than its neighbours (contained dpT reco bin [52.7, 76.0] is
            #      2.91% of the axis against 3.6-8.3% for the rest), and if the cutoff
            #      sits above it that whole row loses its labels and reads as empty.
            #      The default is set below the narrowest bin any current binning
            #      produces; raise it only after checking the reco widths.
            #   2. CELL SIZE: cells narrower/shorter than annotate_min_cell (as a
            #      fraction of the axis span) have no room for text at all. This is what
            #      happens where the binning is necessarily narrow - e.g. the
            #      equal-population bins packed into a sharp forward peak, which can be
            #      1-2% of the axis wide while still holding a large diagonal fraction.
            # In both cases the colour scale still represents the cell faithfully.
            # Set either to 0 to disable that suppression.
            # The visible window may be narrower than the binning: an overflow bin
            # (e.g. truth [p_max, inf)) belongs in the response matrix used by the fit,
            # but plotting it would compress the physical range into a sliver of the
            # axis. Setting xrange/yrange on the artist crops the view to the physical
            # region while the column normalisation - and hence every displayed
            # fraction - still accounts for the overflow.
            vis_x = (float(x_edges[0]), float(x_edges[-1])) if self._xrange is None \
                else (float(self._xrange[0]), float(self._xrange[1]))
            vis_y = (float(y_edges[0]), float(y_edges[-1])) if self._yrange is None \
                else (float(self._yrange[0]), float(self._yrange[1]))
            x_span = vis_x[1] - vis_x[0]
            y_span = vis_y[1] - vis_y[0]
            for i in range(ny):
                for j in range(nx):
                    percentage = normalized_values.T[i, j] * 100
                    if percentage <= annotate_threshold:
                        continue
                    # skip cells lying wholly outside the visible window
                    if x_edges[j + 1] <= vis_x[0] or x_edges[j] >= vis_x[1]:
                        continue
                    if y_edges[i + 1] <= vis_y[0] or y_edges[i] >= vis_y[1]:
                        continue
                    if annotate_min_cell > 0:
                        # clip to the visible window so a bin that merely extends past
                        # the edge is judged on the part actually drawn
                        w = min(x_edges[j + 1], vis_x[1]) - max(x_edges[j], vis_x[0])
                        h = min(y_edges[i + 1], vis_y[1]) - max(y_edges[i], vis_y[0])
                        if w / x_span < annotate_min_cell or h / y_span < annotate_min_cell:
                            continue
                    x_pos = (min(x_edges[j + 1], vis_x[1]) + max(x_edges[j], vis_x[0])) / 2
                    y_pos = (min(y_edges[i + 1], vis_y[1]) + max(y_edges[i], vis_y[0])) / 2
                    text_color = 'white' if percentage > 12.5 else 'black'
                    ax.text(x_pos, y_pos, f'{percentage:.1f}%',
                            ha='center', va='center',
                            color=text_color, fontsize=6, weight='bold')

            ax.set_xlim(*vis_x)
            ax.set_ylim(*vis_y)

            # Set labels and formatting
            ax.set_xlabel(self._variables[0]._xlabel if self._xtitle is None else self._xtitle, fontsize=12,
                          weight='bold')
            ax.set_ylabel(self._variables[1]._xlabel if self._ytitle is None else self._ytitle, fontsize=12,
                          weight='bold')

            ax.tick_params(axis='both', which='major',
                           labelsize=12,
                           size=8,
                           width=2)

            # Automatic locators pick round intervals that need not land on
            # the end of the range (an angular axis on [0, 180] stops
            # labelling at 150), so an explicit xticks list on the Variable
            # overrides them per axis.
            if getattr(self._variables[0], '_xticks', None):
                ax.set_xticks(self._variables[0]._xticks)
            else:
                ax.xaxis.set_major_locator(plt.MaxNLocator(nbins=6, steps=[1,2,5,10]))
            if getattr(self._variables[1], '_xticks', None):
                ax.set_yticks(self._variables[1]._xticks)
            else:
                ax.yaxis.set_major_locator(plt.MaxNLocator(nbins=6, steps=[1,2,5,10]))

            # Add colorbar
            if draw_colorbar:
                cbar = plt.colorbar(im, ax=ax)
                cbar.set_label('Column Normalized Fraction', fontsize=12, weight='bold')
                cbar.ax.tick_params(labelsize=12, size=6, width=1.5)
        if show_option == 'projection' and self._plotdata_diagonal is not None:
            labels, data = zip(*self._plotdata_diagonal.items())
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges_diagonal[l][:-1] + np.diff(self._binedges_diagonal[l]) / 2 for l in labels]

            ax.hist(bincenters, weights=data, bins=self._variables[0]._nbins,
                    range=(-1,1) if self._xrange is None else self._xrange,
                    histtype='barstacked', label=labels, color=colors, stacked=True)
            ax.set_xlabel('(Y-X)/X' if self._xtitle is None else self._xtitle, fontsize=14, weight='bold')
            ax.set_ylabel('Entries', fontsize=12, weight='bold')
            ax.set_xlim(-1, 1) if self._xrange is None else ax.set_xlim(self._xrange)

            # Set tick mark size and tick label font size
            ax.tick_params(axis='both', which='major',
                           labelsize=12,  # Font size of tick labels
                           size=8,  # Length of major tick marks
                           width=2)  # Width/thickness of tick marks

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges_diagonal[labels[0]], fit_type, range=(-1,1) if self._xrange is None else self._xrange)

            if draw_stat_error:
                x = bincenters[0]
                y = np.sum(data, axis=0)
                xerr = np.diff(self._binedges_diagonal[labels[0]]) / 2
                yerr = np.sqrt(y)
                draw_error_boxes(ax, x, y, xerr, yerr, facecolor='gray', edgecolor='none', alpha=0.5, hatch='///')

            # Drop legend entries for categories with no entries in the
            # projection (e.g. the empty 'Data' category on MC-only plots)
            empty_labels = {lab for lab, d in self._plotdata_diagonal.items()
                            if np.sum(d) == 0}

            if invert_stack_order:
                h, l = ax.get_legend_handles_labels()
                filtered = [(h, l) for h, l in zip(h, l)
                            if 'QE' not in l and l not in empty_labels]
                if filtered:
                    h, l = zip(*filtered)
                h, l = list(h), list(l)
                if draw_stat_error:
                    h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                    l.append('MC Statistical Uncertainty')
                ax.legend(h[-2::-1] + h[-1:], l[-2::-1] + l[-1:], fontsize=10)
            else:
                h, l = ax.get_legend_handles_labels()
                filtered = [(hh, ll) for hh, ll in zip(h, l)
                            if 'QE' not in ll and ll not in empty_labels]
                if filtered:
                    h, l = zip(*filtered)
                h, l = list(h), list(l)
                if draw_stat_error:
                    h.append(plt.Rectangle((0, 0), 1, 1, fc='gray', alpha=0.5, hatch='///'))
                    l.append('MC Statistical Uncertainty')
                ax.legend(h, l, fontsize=10)
        
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
        if style.mark_pot and (show_option not in ['smearing', '2d']):
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=0.035 if (style.scilimits and not logy) is not None else 0)

        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            xr = self._variables[0]._range if self._xrange is None else self._xrange
            if xr[0] == 0:    
                xhigh_exporder = np.floor(np.log10(xr[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, xr[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')