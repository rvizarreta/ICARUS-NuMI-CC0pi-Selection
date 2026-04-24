import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from artists import SpineArtist
from style import Style
from variable import Variable
from utilities import mark_pot, mark_preliminary


class SpinePurity(SpineArtist):
    """
    SpinePurity class optimized for showing cut progression using signal tree.
    """

    def __init__(self, variable, categories, cuts, title,
                 xrange=None, xtitle=None, show_option='table',
                 signal_categories=[0.0], npts=1e6):
        print(f"SpinePurity.__init__: Creating signal tree purity artist")
        print(f"   - Variable: {variable._key}")
        print(f"   - Signal categories: {signal_categories}")

        super().__init__()
        self._variable = variable
        self._categories = categories
        self._title = title
        self._cuts = cuts
        self._xrange = xrange
        self._xtitle = xtitle
        self._show_option = show_option
        self._signal_categories = signal_categories
        self._npts = int(npts)
        self._purity_results = {}
        self._samples = []

    def add_sample(self, sample, is_ordinate=False):
        print(f"SpinePurity.add_sample: Adding sample {sample._name}")
        self._samples.append(sample)
        self.calculate_progression(sample)

    def calculate_progression(self, sample):
        """Calculate purity progression using both signal and selected trees."""
        print(f"SpinePurity.calculate_progression: Starting calculation")

        # Access the file handle to get both trees
        file_handle = sample._file_handle

        try:
            # Load both trees
            signal_tree = file_handle['signal']
            selected_tree = file_handle['selected']

            signal_df = signal_tree.arrays(library='pd')
            selected_df = selected_tree.arrays(library='pd')

            print(f"   - Signal tree events: {len(signal_df)}")
            print(f"   - Selected tree events: {len(selected_df)}")

            # Calculate purity progression
            stages = {}

            # Stage 1: Signal tree (initial purity)
            signal_events_in_signal = len(signal_df[signal_df['true_category'].isin(self._signal_categories)])
            total_signal_events = len(signal_df)
            signal_purity = signal_events_in_signal / total_signal_events if total_signal_events > 0 else 0.0
            stages['initial_selection'] = signal_purity

            print(f"   1. Initial Selection (Signal Tree): {signal_purity:.3f} ({signal_purity * 100:.1f}%)")
            print(f"      Total: {total_signal_events}, Signal: {signal_events_in_signal}")

            # Show intermediate stages if we can simulate them from signal tree
            intermediate_stages = self.simulate_cut_progression(signal_df)
            stages.update(intermediate_stages)

            # Final stage: Selected tree (final purity)
            signal_events_in_selected = len(selected_df[selected_df['true_category'].isin(self._signal_categories)])
            total_selected_events = len(selected_df)
            selected_purity = signal_events_in_selected / total_selected_events if total_selected_events > 0 else 0.0
            stages['final_selection'] = selected_purity

            print(f"   N. Final Selection (Selected Tree): {selected_purity:.3f} ({selected_purity * 100:.1f}%)")
            print(f"      Total: {total_selected_events}, Signal: {signal_events_in_selected}")

            # Calculate efficiency
            efficiency = total_selected_events / total_signal_events if total_signal_events > 0 else 0.0
            print(f"   Selection Efficiency: {efficiency:.3f} ({efficiency * 100:.1f}%)")

            # Store results
            self._purity_results["Selection"] = stages

        except Exception as e:
            print(f"   Error accessing trees: {e}")
            # Fallback to single tree
            self.calculate_single_tree(sample)

    def simulate_cut_progression(self, signal_df):
        """Simulate cut progression using signal tree data."""
        print(f"   - Simulating cut progression from signal tree...")

        stages = {}
        current_df = signal_df.copy()
        stage_num = 2

        # Check what columns are available for cuts
        available_cuts = []

        # Define logical cut progression based on typical 1μNp0π selection
        potential_cuts = [
            ('fiducial', 'reco_fiducial == 1', 'Fiducial Volume'),
            ('containment', 'reco_containment_cut == 1', 'Containment'),
            ('single_muon', 'reco_muon_multiplicity == 1', 'Single Muon'),
            ('proton_mult', 'reco_proton_multiplicity >= 1', '≥1 Proton'),
            ('no_pions', 'reco_pion_multiplicity == 0', 'No Charged Pions'),
            ('no_photons', 'reco_photon_multiplicity == 0', 'No Photons'),
            ('no_electrons', 'reco_electron_multiplicity == 0', 'No Electrons')
        ]

        for cut_name, cut_expression, cut_label in potential_cuts:
            try:
                # Test if we can apply this cut
                cut_mask = current_df.eval(cut_expression)
                current_df = current_df[cut_mask]

                # Calculate purity after this cut
                signal_events = len(current_df[current_df['true_category'].isin(self._signal_categories)])
                total_events = len(current_df)
                purity = signal_events / total_events if total_events > 0 else 0.0

                stages[cut_name] = purity
                print(f"   {stage_num}. {cut_label}: {purity:.3f} ({purity * 100:.1f}%) - {total_events} events")
                stage_num += 1

            except Exception as e:
                print(f"      Warning: Could not apply {cut_label}: {e}")
                continue

        return stages

    def calculate_single_tree(self, sample):
        """Fallback single tree calculation."""
        print(f"   - Fallback to single tree calculation")

        data = sample._data
        signal_events = len(data[data['true_category'].isin(self._signal_categories)])
        total_events = len(data)
        purity = signal_events / total_events if total_events > 0 else 0.0

        self._purity_results["Selection"] = {'final_selection': purity}

    def reduce(self, group, significance=0.6827):
        print(f"SpinePurity.reduce: Called for group '{group}'")

        if "Selection" not in self._purity_results:
            raise ValueError("No purity results available")

        purity_dict = self._purity_results["Selection"]

        # Convert to efficiency-like format
        seq_results = {}
        unseq_results = {}

        for i, (cut_key, purity_value) in enumerate(purity_dict.items()):
            seq_key = f'unbinned_seq_{cut_key}'
            unseq_key = f'unbinned_unseq_{cut_key}'
            seq_results[seq_key] = purity_value
            unseq_results[unseq_key] = purity_value

        combined_results = {**seq_results, **unseq_results}
        return None, combined_results, {}, {}

    def draw(self, ax, show_option='table', percentage=True, groups=None, **kwargs):
        print(f"SpinePurity.draw: Creating progression table")

        if show_option != 'table':
            raise NotImplementedError("Only 'table' option is currently supported")

        if "Selection" not in self._purity_results:
            raise ValueError("No purity results available")

        # Clear axis
        ax.axis('off')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # Add title
        title = self._title or "1μNp0π Selection Purity Progression"
        ax.set_title(title, pad=20, fontsize=16, fontweight='bold')

        # Format data
        formatter = (lambda x: f'{100 * x:.1f}%') if percentage else (lambda x: f'{x:.3f}')
        purity_key = 'Purity [%]' if percentage else 'Purity'

        # Get purity results
        purity_dict = self._purity_results["Selection"]

        # Create meaningful labels
        stage_labels = {
            'initial_selection': 'Initial Selection (Signal Tree)',
            'fiducial': '+ Fiducial Volume',
            'containment': '+ Containment',
            'single_muon': '+ Single Muon',
            'proton_mult': '+ ≥1 Proton',
            'no_pions': '+ No Charged Pions',
            'no_photons': '+ No Photons',
            'no_electrons': '+ No Electrons',
            'final_selection': 'Final Selection (Selected Tree)'
        }

        # Create table data
        cuts_list = []
        purities_list = []

        # Order stages logically
        stage_order = ['initial_selection', 'fiducial', 'containment', 'single_muon',
                       'proton_mult', 'no_pions', 'no_photons', 'no_electrons', 'final_selection']

        for stage in stage_order:
            if stage in purity_dict:
                purity_value = purity_dict[stage]
                label = stage_labels.get(stage, stage.replace('_', ' ').title())

                cuts_list.append(label)
                purities_list.append(formatter(purity_value))

        # Create table
        table_data = list(zip(cuts_list, purities_list))
        headers = ['Selection Stage', purity_key]

        print(f"   - Creating table with {len(table_data)} rows")

        table = ax.table(cellText=table_data, colLabels=headers,
                         loc='center', cellLoc='center',
                         bbox=[0.05, 0.2, 0.9, 0.6])

        # Style the table
        table.scale(1, 2.0)
        table.set_fontsize(11)

        # Color and style cells
        for (i, j), cell in table.get_celld().items():
            if i == 0:  # Header row
                cell.set_text_props(weight='bold', fontsize=13)
                cell.set_facecolor('#1976D2')  # Blue header
                cell.set_text_props(color='white')
            else:
                # Light blue for data rows
                cell.set_facecolor('#E3F2FD')

        # Add summary
        if 'final_selection' in purity_dict:
            final_purity = purity_dict['final_selection']
            summary_text = f"Final Purity: {formatter(final_purity)}\nProgression: Signal Tree → Selected Tree"
            ax.text(0.5, 0.1, summary_text, ha='center', va='center',
                    transform=ax.transAxes, fontsize=12, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.5", facecolor="#4CAF50", alpha=0.9, edgecolor='black'))

        return table

    def get_purity_results(self):
        """Get the calculated purity results."""
        return self._purity_results