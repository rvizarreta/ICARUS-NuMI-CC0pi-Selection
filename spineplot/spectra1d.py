import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from spectra import SpineSpectra
from style import Style
from variable import Variable
from systematic import Systematic
from utilities import mark_pot, mark_preliminary, draw_error_boxes

class SpineSpectra1D(SpineSpectra):
    """
    A class designed to encapsulate a single variable's spectrum for an
    ensemble of samples. This is a specialization of the SpineSpectra
    class that is intended to be used for 1D spectra. At its core, this
    is a simple histogram plot of a single variable.

    Attributes
    ----------
    _title : str
        The title of the spectrum. This will be placed at the top of
        the axis assigned to the artist.
    _xrange : tuple
        The range of the x-axis for the spectrum. This is a tuple of
        the lower and upper limits of the x-axis. If None, the range
        will be determined by the range set in the Variable object.
    _xtitle : str
        The label for the x-axis. If None, the label will be taken
        from the Variable object.
    _yrange : tuple, or float, optional
        If this is a tuple, it is the range of the y-axis for the
        spectrum. If this is a float, it will scale the maximum value
        of the histogram by this factor. If None, the range will be
        determined by the range of the histogram.
    _ytitle : str
        The label for the y-axis. If None, the label will be 'Candidates'.
    _variable : Variable
        The Variable object for the spectrum.
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
    def __init__(self, variable, categories, colors, category_types,
                 title=None, xrange=None, xtitle=None,
                 yrange=None, ytitle=None) -> None:
        """
        Initializes the SpineSpectra1D.

        Parameters
        ----------
        variable : Variable
            The Variable object for the spectrum.
        categories : dict
            A dictionary of the categories for the spectrum. This
            serves as a map between the category label in the input
            TTree and the category label for the spectrum (and
            therefore what is shown in a single legend entry).
        colors : dict
            A dictionary of the colors for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the color to use for the histogram. The color can be
            any valid matplotlib color string or a cycle indicator
            (e.g. 'C0', 'C1', etc.).
        category_types : dict
            A dictionary of the types for the categories in the
            spectrum. This serves as a map between the category label
            for the spectrum (value in the `_categories` dictionary)
            and the type of plot to use for the histogram. The type
            should be either 'histogram' or 'scatter' to correspond to
            a stacked histogram or scatter plot, respectively.
        title : str, optional
            The title of the spectrum. This will be placed at the top
            of the axis assigned to the artist. The default is None.
        xrange : tuple, optional
            The range of the x-axis for the spectrum. This is a tuple
            of the lower and upper limits of the x-axis. If None, the
            range will be determined by the range set in the Variable
            object. The default is None.
        xtitle : str, optional
            The label for the x-axis. If None, the label will be taken
            from the Variable object. The default is None.
        yrange : tuple, or float, optional
            If this is a tuple, it is the range of the y-axis for the
            spectrum. If this is a float, it will scale the maximum
            value of the histogram by this factor. If None, the range
            will be determined by the range of the histogram. The
            default is None.
        ytitle : str, optional
            The label for the y-axis. If None, the label will be
            'Candidates'. The default is None.

        Returns
        -------
        None.
        """
        super().__init__([variable,], categories, colors, title,
                         xrange, xtitle, yrange, ytitle)
        self._variable = self._variables[0]
        self._category_types = category_types

    def add_sample(self, sample, is_ordinate) -> None:
        """
        Adds a sample to the SpineSpectra1D object. The sample's data
        is extracted per category and stored for later plotting.
        Multiple samples may have overlapping categories, so the data
        is stored in a dictionary with the category as the key.

        Parameters
        ----------
        sample : Sample
            The sample to add to the SpineSpectra1D object.
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
            self._onebincount = {}
        data, weights = sample.get_data([self._variable._key,], with_mask=self._variable.mask)
        for category, values in data.items():
            values = values[0]
            if category not in self._categories.keys():
                continue
            if self._categories[category] not in self._plotdata:
                self._plotdata[self._categories[category]] = np.zeros(self._variable._nbins)
                self._onebincount[self._categories[category]] = 0
            xr = self._variable._range if self._xrange is None else self._xrange
            h = np.histogram(values, bins=self._variable._nbins, range=xr, weights=weights[category])
            self._onebincount[self._categories[category]] += np.sum(weights[category])
            self._plotdata[self._categories[category]] += h[0]
            self._binedges[self._categories[category]] = h[1]

    def draw(self, ax, style, show_component_number=False,
             show_component_percentage=False, invert_stack_order=False,
             fit_type=None, logx=False, logy=False, normalize=False,
             draw_error=None, draw_ratio=False, prediction=True, show_fraction=False) -> None:
        """
        Plots the data for the SpineSpectra1D object.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the artist on.
        style : Style
            The style to use when drawing the artist. The default is
            None. This is intended to be used in cases where the artist
            has some configurable style options.
        show_component_number : bool
            A flag to indicate if the component number should be shown
            in the legend. The default is False.
        show_component_percentage : bool
            A flag to indicate if the component percentage should be
            shown in the legend. The default is False.
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
        normalize : bool
            A flag to indicate if the histogram should be normalized
            to unity. The default is False.
        draw_error : str, optional
            Indicates the name of the Systematic object to use for
            drawing the error boxes. The default is None.
        prediction : bool, optional
            A flag to indicate if the prediction line legend (sum of all MC
            components) should be drawn. The default is True.
        show_fraction : bool, optional
            If True, shows the per-bin fraction of each component
            (each bin sums to 1.0). The y-axis becomes 'Fraction'
            instead of 'Candidates'. The default is False.

        Returns
        -------
        None.
        """
        # Handle ratio plot setup
        if draw_ratio:
            # Split the current axes into two subplots
            fig = ax.figure
            gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.0)
            ax.remove()  # Remove the original axis
            ax_main = fig.add_subplot(gs[0])
            ax_ratio = fig.add_subplot(gs[1], sharex=ax_main)
            ax = ax_main  # Use ax_main for the rest of the drawing
        else:
            ax_ratio = None

        ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle, fontsize=12, weight='bold')
        ylabel = 'Fraction' if show_fraction else 'Candidates'
        ax.set_ylabel(ylabel, fontsize=12, weight='bold')
        ax.set_xlim(*self._variable._range if self._xrange is None else self._xrange)
        ax.set_title(self._title)

        # Set tick mark size and tick label font size
        ax.tick_params(axis='both', which='major',
                       labelsize=12,  # Font size of tick labels
                       size=8,  # Length of major tick marks
                       width=2)  # Width/thickness of tick marks

        if self._plotdata is not None:
            labels, data = zip(*self._plotdata.items())
            # Store original labels before modification
            original_labels = list(labels)
            colors = [self._colors[label] for label in labels]
            bincenters = [self._binedges[l][:-1] + np.diff(self._binedges[l]) / 2 for l in labels]
            binwidths = [np.diff(self._binedges[l]) for l in labels]
            xr = self._variable._range if self._xrange is None else self._xrange

            histogram_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'histogram']
            scatter_mask = [li for li, label in enumerate(labels) if self._category_types[label] == 'scatter']

            denominator = np.sum([self._onebincount[labels[i]] for i in histogram_mask])
            counts = [x for x in self._onebincount.values()]

            if fit_type is not None:
                super().fit_with_function(ax, bincenters[0], np.sum(data, axis=0), self._binedges[labels[0]], fit_type, range=xr)

            if show_component_number and show_component_percentage:
                # Find Signal and Signal QE indices
                signal_idx = None
                signal_qe_idx = None
                for li, label in enumerate(labels):
                    if label == 'Signal':
                        signal_idx = li
                    elif label == 'Signal QE':
                        signal_qe_idx = li

                # Calculate combined signal
                if signal_idx is not None and signal_qe_idx is not None:
                    signal_events = np.sum(counts[signal_idx])
                    signal_qe_events = np.sum(counts[signal_qe_idx])
                    combined_signal = signal_events + signal_qe_events

                # Build labels with custom logic
                new_labels = []
                for li, (label, d) in enumerate(zip(labels, counts)):
                    if li in histogram_mask:
                        if label == 'Signal' and signal_idx is not None and signal_qe_idx is not None:
                            # Signal: show combined total
                            new_labels.append(f'{label} ({combined_signal:.1f}, {combined_signal / denominator:.2%})')
                        elif label == 'Signal QE' and signal_idx is not None and signal_qe_idx is not None:
                            # Signal QE: show what fraction of combined signal is QE
                            new_labels.append(
                                f'{label} ({signal_qe_events / combined_signal:.2%} of signal)')
                        else:
                            # All other histogram components
                            new_labels.append(f'{label} ({np.sum(d):.1f}, {np.sum(d) / denominator:.2%})')
                    else:
                        # Scatter components
                        new_labels.append(f'{label} ({np.sum(d):.1f})')

                labels = new_labels

            elif show_component_number:
                labels = [f'{label} ({np.sum(d):.1f})' for label, d in zip(labels, counts)]
            elif show_component_percentage:
                labels = [f'{label} ({np.sum(d) / denominator:.2%})' if li in histogram_mask else label for
                          li, (label, d) in enumerate(zip(labels, counts))]

            if invert_stack_order:
                reduce = lambda x : [x[i] for i in histogram_mask[::-1]]
            else:
                reduce = lambda x : [x[i] for i in histogram_mask]

            scale = 1.0 if not normalize else 1.0 / np.sum(reduce(data))

            # SAVE MC prediction for ratio plot BEFORE any modifications
            mc_sum_for_ratio = np.sum(reduce(data), axis=0)

            # Apply per-bin normalization AFTER reduce() if show_fraction is True
            if show_fraction:
                # Get only the histogram components
                hist_data = reduce(data)
                # Get the total per bin for histogram components only
                total_per_bin = np.sum(hist_data, axis=0)
                # Avoid division by zero
                total_per_bin = np.where(total_per_bin > 0, total_per_bin, 1.0)
                # Normalize each histogram component by the total in each bin
                hist_data_normalized = [d / total_per_bin for d in hist_data]
                # We need to modify the reduce function to return normalized data
                # Store normalized histogram data
                normalized_hist_dict = {labels[histogram_mask[i]]: hist_data_normalized[i]
                                        for i in range(len(histogram_mask))}
                # Create new data tuple with normalized histogram data
                data_list = list(data)
                for idx, i in enumerate(histogram_mask):
                    data_list[i] = hist_data_normalized[
                        idx if not invert_stack_order else len(histogram_mask) - 1 - idx]
                data = tuple(data_list)
                scale = 1.0

            if draw_error:
                systs = [s[draw_error] for s in self._systematics.values() if draw_error in s]
                cov = np.sum(s.get_covariance(self._variable._key) for s in systs)
                x = reduce(bincenters)[0]
                y = scale * np.sum(reduce(data), axis=0)
                xerr = [x / 2 for x in binwidths[0]]
                #scov = Systematic.transform_as(cov, scale if not normalize else np.sum(reduce(data), axis=0))
                scov = Systematic.transform_as(cov, scale)
                yerr = np.sqrt(np.diag(scov))

                draw_error_boxes(ax, x, y, xerr, yerr, facecolor='lightgray', edgecolor='gray', alpha=1.0, hatch='xxx', linewidth=0.0)

                # Draw border line on top of uncertainty band
                bin_edges = list(self._binedges.values())[0]
                y_top = y + yerr
                ax.step(bin_edges,
                        np.append(y_top, y_top[-1]),
                        where='post',
                        color='grey',  # or 'black', 'gray', etc.
                        linewidth=1.5,
                        zorder=1)

            if show_fraction:
                # Draw each component as a separate line (not stacked)
                for i, (label, d, color, bc) in enumerate(
                        zip(reduce(labels), reduce(data), reduce(colors), reduce(bincenters))):
                    bin_edges = self._binedges[
                        original_labels[histogram_mask[i] if not invert_stack_order else histogram_mask[::-1][i]]]
                    ax.step(bin_edges, np.append(scale * d, scale * d[-1]),
                            where='post', label=label, color=color, linewidth=1.5)
            else:
                # Draw stacked histogram (original behavior)
                ax.hist(reduce(bincenters), weights=[scale * x for x in reduce(data)], bins=self._variable._nbins,
                        range=xr, label=reduce(labels), color=reduce(colors), alpha=0.85, **style.plot_kwargs)

            # Draw top border line for components with 'QE' in label (only for stacked plots)
            if not show_fraction:
                cumulative_heights = np.zeros(self._variable._nbins)
                for i, (label, d, orig_label) in enumerate(zip(reduce(labels), reduce(data), reduce(original_labels))):
                    cumulative_heights += scale * d

                    # Check if this component has 'QE' in label
                    if 'QE' in label:
                        # Use the ORIGINAL label to access binedges
                        bin_edges = self._binedges[orig_label]  # <-- Use orig_label instead
                        ax.step(bin_edges,
                                np.append(cumulative_heights, cumulative_heights[-1]),
                                where='post',
                                color='darkred',
                                linewidth=2,
                                zorder=10)
                        break

            # Add prediction line (sum of all MC components) - skip for fraction plots
            if not show_fraction:
                total_prediction = scale * np.sum(reduce(data), axis=0)
                total_events = np.sum(reduce(data))  # Sum all events in histogram categories
                first_hist_category = list(self._plotdata.keys())[histogram_mask[0]]
                bin_edges = self._binedges[first_hist_category]
                label = f'Prediction ({total_events:.1f} events)' if prediction else None
                ax.step(bin_edges, np.append(total_prediction, total_prediction[-1]),
                        where='post', color='black', linewidth=1.5, label=label)


            reduce = lambda x : [x[i] for i in scatter_mask]
            for i, label in enumerate(reduce(labels)):
                idx = scatter_mask[i]  # Get the actual index in the original arrays
                scale = 1.0 if not normalize else 1.0 / np.sum(data[idx])

                # Filter out bins with no data
                data_mask = data[scatter_mask[i]] > 0
                if not np.any(data_mask):  # Skip if no data at all
                    continue

                xerr_data = [x / 2 for x in binwidths[idx]]  # Get bin widths for this category
                ax.errorbar(bincenters[scatter_mask[i]][data_mask],
                            scale * data[scatter_mask[i]][data_mask],
                            xerr=[xerr_data[j] for j in range(len(xerr_data)) if data_mask[j]],
                            yerr=scale * np.sqrt(data[scatter_mask[i]][data_mask]),
                            fmt='o',
                            markersize=6,
                            markerfacecolor='black',
                            markeredgecolor='black',
                            color='black',
                            capsize=3,
                            capthick=1.5,
                            elinewidth=1.5,
                            label=label)
        if invert_stack_order:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='lightgray', ec='gray', alpha=0.6, hatch='xxx', linewidth=0.8))
                l.append(systs[0].label)
                ax.legend(h[-2::-1]+h[-1:], l[-2::-1]+l[-1:], fontsize=8)
            else:
                ax.legend(h[::-1], l[::-1], fontsize=8)
        else:
            h, l = ax.get_legend_handles_labels()
            if draw_error:
                h.append(plt.Rectangle((0, 0), 1, 1, fc='lightgray', ec='gray', alpha=0.6, hatch='xxx', linewidth=0.8))
                l.append(systs[0].label)
            ax.legend(h, l, fontsize=8)

        # Add borders to legend patches for 'QE' entries
        legend = ax.get_legend()
        for handle, label_text in zip(legend.legend_handles, legend.get_texts()):
            if 'QE' in label_text.get_text():
                handle.set_edgecolor('darkred')
                handle.set_linewidth(1)
                handle.set_facecolor('white')

        # Automatically extend y-axis by 35% to make room for legend
        #if self._yrange is None:
        #    current_ylim = ax.get_ylim()
        #    ax.set_ylim(current_ylim[0], current_ylim[1] * 1.35)

        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            yl = ax.get_ylim()[1]
            ax.set_ylim(None, yl * self._yrange)

        if draw_ratio and ax_ratio is not None and scatter_mask:
            # Get data from self._plotdata directly
            data_values = self._plotdata['Data']
            data_centers = bincenters[scatter_mask[0]]

            # Use saved MC prediction
            mc_prediction_ratio = mc_sum_for_ratio
            mc_stat_err = np.sqrt(mc_prediction_ratio)

            # Calculate data/mc ratio
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.where(data_values > 0, data_values / mc_prediction_ratio, 0)

            # Calculate uncertainties
            mc_err = np.zeros_like(mc_prediction_ratio)
            if draw_error:
                systs_ratio = [s[draw_error] for s in self._systematics.values() if draw_error in s]
                if systs_ratio:
                    cov_ratio = np.sum(s.get_covariance(self._variable._key) for s in systs_ratio)
                    scov_ratio = Systematic.transform_as(cov_ratio, 1.0)
                    mc_syst_err = np.sqrt(np.diag(scov_ratio))
                    # Add MC statistical and systematic errors in quadrature
                    mc_err = np.sqrt(mc_stat_err ** 2 + mc_syst_err ** 2)
                else:
                    # If no systematics, just use MC statistical error
                    mc_err = mc_stat_err
            else:
                # If draw_error is False, still include MC statistical error in band
                mc_err = mc_stat_err

            data_err = np.sqrt(data_values)

            with np.errstate(divide='ignore', invalid='ignore'):
                ratio_err = np.where(mc_prediction_ratio > 0,
                                     data_err / mc_prediction_ratio, 0)

            # Draw ratio plot
            ax_ratio.axhline(y=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)

            if draw_error and np.any(mc_err > 0):
                with np.errstate(divide='ignore', invalid='ignore'):
                    rel_mc_err = np.where(mc_prediction_ratio > 0, mc_err / mc_prediction_ratio, 0)
                # Only draw error boxes where data exists
                if np.any(data_mask):
                    xerr_ratio = np.array(binwidths[scatter_mask[0]])[data_mask] / 2
                    draw_error_boxes(ax_ratio, data_centers[data_mask], np.ones(np.sum(data_mask)),
                                     xerr_ratio, rel_mc_err[data_mask],
                                     facecolor='lightgray', edgecolor='gray', alpha=0.4,
                                     hatch='xxx', linewidth=0.6)

            # Only plot data points where data exists
            if np.any(data_mask):
                xerr_ratio_filtered = np.array(binwidths[scatter_mask[0]])[data_mask] / 2
                ax_ratio.errorbar(data_centers[data_mask], ratio[data_mask],
                                  xerr=xerr_ratio_filtered,
                                  yerr=ratio_err[data_mask],
                                  fmt='o', markersize=4,
                                  markerfacecolor='black', markeredgecolor='black',
                                  color='black', capsize=2, elinewidth=1)

            ax_ratio.set_ylabel('Data/MC', fontsize=10, weight='bold')
            ax_ratio.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle,
                                fontsize=12, weight='bold')
            ax_ratio.set_ylim(0.0, 2.0)
            ax_ratio.grid(True, alpha=0.3, axis='y')
            ax_ratio.tick_params(axis='both', which='major', labelsize=10)
            ax_ratio.set_xlim(*self._variable._range if self._xrange is None else self._xrange)

            ax.set_xlabel('')
            ax.tick_params(labelbottom=False)

        # Set the axis to be logarithmic if requested.
        if logx:
            # Modify the x-axis limits to ensure that the lower limit
            # is greater than zero. The lower edge needs to be at least
            # 3 orders of magnitude less than the maximum value in the
            # plot.
            xr = self._variable._range if self._xrange is None else self._xrange
            if xr[0] == 0:    
                xhigh_exporder = np.floor(np.log10(xr[1]))
                xlow = xhigh_exporder - 3
                ax.set_xlim(10**xlow, xr[1])
            ax.set_xscale('log')
        if logy:
            ax.set_yscale('log')

        # hadj and vadj are used to adjust the position of the POT and
        # preliminary labels horizontally and vertically, respectively.
        # This is necessary to ensure that the labels do not overlap
        # with plot elements. The following logic is meant to capture
        # all cases where the labels might overlap with the plot.
        hadj = 0
        vadj = 0

        # The scilimits option cannot be used with a logarithmic y-axis.
        # The hadj value is adjusted to ensure that the POT label does
        # not overlap with the scientific notation placed above the
        # y-axis.
        if style.scilimits and not logy:
            ax.ticklabel_format(axis='y', scilimits=style.scilimits)
            hadj = 0.035

        # The vadj value is adjusted to ensure that the POT label does
        # not overlap with the top axis of the plot when the y-axis is
        # logarithmic.
        if logy:
            vadj = 0.1
        
        # Add the POT and preliminary labels to the plot.
        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal, vadj=vadj)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary, hadj=hadj, vadj=vadj)

    def draw_systematics(self, ax, style, systematic_names=None,
                         fractional=True, show_total=True,
                         show_individual=True) -> None:
        """
        Plots systematic uncertainties as a function of the variable.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to draw the systematics on.
        style : Style
            The style to use when drawing the plot.
        systematic_names : list of str, optional
            List of systematic names to plot. If None, plots all available
            systematics. The default is None.
        fractional : bool, optional
            If True, plot fractional uncertainties (σ/y). If False, plot
            absolute uncertainties (σ). The default is True.
        show_total : bool, optional
            If True, show the total combined systematic uncertainty.
            The default is True.
        show_individual : bool, optional
            If True, show individual systematic uncertainties. The default
            is True.
        """
        # Set axis labels
        ax.set_xlabel(self._variable._xlabel if self._xtitle is None else self._xtitle, fontsize=12, weight='bold')
        if fractional:
            ax.set_ylabel('Fractional Uncertainty', fontsize=12, weight='bold')
        else:
            ax.set_ylabel('Absolute Uncertainty', fontsize=12, weight='bold')
        ax.set_xlim(*self._variable._range if self._xrange is None else self._xrange)

        # Get bin centers for plotting
        if self._plotdata is None:
            print("Warning: No data has been added to the spectrum yet.")
            return

        labels = list(self._plotdata.keys())
        bincenters = self._binedges[labels[0]][:-1] + np.diff(self._binedges[labels[0]]) / 2

        # Get the prediction (central value) for normalization
        histogram_mask = [li for li, label in enumerate(labels)
                          if self._category_types[label] == 'histogram']
        data = list(self._plotdata.values())
        y_pred = np.sum([data[i] for i in histogram_mask], axis=0)

        # Determine which systematics to plot
        if systematic_names is None:
            # Get all systematics from all samples
            all_syst_names = set()
            for sample_systs in self._systematics.values():
                all_syst_names.update(sample_systs.keys())
            systematic_names = sorted(list(all_syst_names))

        # Define colors for individual systematics
        #colors = ['green', 'mediumblue', 'darkred', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
        colors = [
            'blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'orange', 'purple', 'brown',
            'aliceblue', 'antiquewhite', 'aqua', 'aquamarine', 'azure', 'beige', 'bisque', 'blanchedalmond',
            'blueviolet', 'burlywood',
            'cadetblue', 'chartreuse', 'chocolate', 'coral', 'cornflowerblue', 'cornsilk', 'crimson', 'darkblue',
            'darkcyan', 'darkgoldenrod',
            'darkgray', 'darkgreen', 'darkkhaki', 'darkmagenta', 'darkolivegreen', 'darkorange', 'darkorchid',
            'darkred', 'darksalmon', 'darkseagreen',
            'darkslateblue', 'darkslategray', 'darkturquoise', 'darkviolet', 'deeppink', 'deepskyblue', 'dimgray',
            'dodgerblue', 'firebrick', 'floralwhite',
            'forestgreen', 'fuchsia', 'gainsboro', 'ghostwhite', 'gold', 'goldenrod', 'gray', 'greenyellow', 'honeydew',
            'hotpink',
            'indianred', 'indigo', 'ivory', 'khaki', 'lavender', 'lavenderblush', 'lawngreen', 'lemonchiffon',
            'lightblue', 'lightcoral',
            'lightcyan', 'lightgoldenrodyellow', 'lightgray', 'lightgreen', 'lightpink', 'lightsalmon', 'lightseagreen',
            'lightskyblue', 'lightslategray', 'lightsteelblue',
            'lightyellow', 'lime', 'limegreen', 'linen', 'maroon', 'mediumaquamarine', 'mediumblue', 'mediumorchid',
            'mediumpurple', 'mediumseagreen',
            'mediumslateblue', 'mediumspringgreen', 'mediumturquoise', 'mediumvioletred', 'midnightblue', 'mintcream',
            'mistyrose', 'moccasin', 'navajowhite', 'navy'
        ]

        # Storage for total uncertainty calculation
        total_cov = None
        stat_cov = None

        # Plot individual systematics
        if show_individual:
            for idx, syst_name in enumerate(systematic_names):
                # Collect systematics from all samples with this name
                systs = []
                for sample_systs in self._systematics.values():
                    if syst_name in sample_systs:
                        systs.append(sample_systs[syst_name])

                if not systs:
                    continue

                # Sum covariance matrices from all samples
                cov = np.sum([s.get_covariance(self._variable._key) for s in systs], axis=0)

                # Check if this is statistical uncertainty
                is_stat = 'stat' in syst_name.lower()

                # Track statistical uncertainty separately
                if is_stat:
                    if stat_cov is None:
                        stat_cov = cov.copy()
                    else:
                        stat_cov += cov

                # Get uncertainty from diagonal
                uncertainty = np.sqrt(np.diag(cov))
                if fractional:
                    print(
                        f"Fractional - Min: {(uncertainty / y_pred).min():.2%}, Max: {(uncertainty / y_pred).max():.2%}")

                # Convert to fractional if requested
                if fractional:
                    # Avoid division by zero
                    with np.errstate(divide='ignore', invalid='ignore'):
                        uncertainty = np.where(y_pred > 0, uncertainty / y_pred, 0)

                # Plot - use dots for statistical, lines for others
                #label = systs[0].label if systs[0].label is not None else syst_name
                base_label = systs[0].label if systs[0].label is not None else syst_name
                mean_frac = np.mean(uncertainty) if fractional else np.mean(
                    np.where(y_pred > 0, uncertainty / y_pred, 0))
                label = f'{base_label} ({mean_frac:.1%})'

                if is_stat:
                    # Plot statistical as black dots
                    ax.step(bincenters, uncertainty, where='mid', label=label, linestyle=':',
                            color='black', linewidth=1)
                else:
                    # Plot systematic as colored line
                    ax.step(bincenters, uncertainty, where='mid', label=label,
                            color=colors[idx], linewidth=1.5, alpha=0.7)

                # Add to total
                if total_cov is None:
                    total_cov = cov.copy()
                else:
                    total_cov += cov

        # Plot total systematic (dashed line for syst only, solid line for total with stat)
        if show_total and total_cov is not None:
            total_uncertainty = np.sqrt(np.diag(total_cov))

            if fractional:
                with np.errstate(divide='ignore', invalid='ignore'):
                    total_uncertainty = np.where(y_pred > 0, total_uncertainty / y_pred, 0)

            # Determine if we have statistical component
            if stat_cov is not None:

                # Plot systematic-only as dashed line
                syst_only_cov = total_cov - stat_cov
                syst_only_uncertainty = np.sqrt(np.diag(syst_only_cov))
                if fractional:
                    with np.errstate(divide='ignore', invalid='ignore'):
                        syst_only_uncertainty = np.where(y_pred > 0, syst_only_uncertainty / y_pred, 0)

                syst_mean = np.mean(syst_only_uncertainty)
                ax.step(bincenters, syst_only_uncertainty, where='mid', label=f'Syst. ({syst_mean:.1%})',
                        color='black', linewidth=1.5, linestyle='--')

                # Total includes stat, so plot as solid line
                total_mean = np.mean(total_uncertainty)
                ax.step(bincenters, total_uncertainty, where='mid', label=f'Total Uncertainty ({total_mean:.1%})',
                        color='black', linewidth=1.5, linestyle='-')
            else:
                # No stat component, just plot total as dashed
                total_mean = np.mean(total_uncertainty)
                ax.step(bincenters, total_uncertainty, where='mid', label=f'Total Systematic ({total_mean:.1%})',
                        color='black', linewidth=2.5, linestyle='--')

        # Add legend
        ax.legend(loc='best', fontsize=10)
        # Add legend outside the plot to the right
        # ax.legend(
        #     loc='upper left',
        #     bbox_to_anchor=(1.02, 1.0),
        #     fontsize=8,
        #     borderaxespad=0.,
        #     ncol=2  # Optional: Split the legend into 2 columns if it is too tall
        # )

        # Add grid for readability
        ax.grid(True, alpha=0.3, linestyle='--')

        # Set y-axis limits if provided
        if isinstance(self._yrange, (tuple, list)):
            ax.set_ylim(*self._yrange)
        elif isinstance(self._yrange, (int, float)):
            yl = ax.get_ylim()[1]
            ax.set_ylim(None, yl * self._yrange)

        # Add POT and preliminary labels if needed
        if style.mark_pot:
            mark_pot(ax, self._exposure, style.mark_pot_horizontal)
        if style.mark_preliminary is not None:
            mark_preliminary(ax, style.mark_preliminary)


