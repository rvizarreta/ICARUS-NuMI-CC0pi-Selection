import os
import sys

workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
sys.path.insert(0, workspace_root)

from analysis import Analysis

DATA_FILE = '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root'

BASE = os.path.dirname(os.path.abspath(__file__))

CONFIGS = {
    'main':     f'{BASE}/leading_proton_momentum.toml',
    'xsec':     f'{BASE}/leading_proton_momentum_xSecUncertainty.toml',
    'xsec_exp': f'{BASE}/leading_proton_momentum_xSecUncertainty_Exp.toml',
    'flux':     f'{BASE}/leading_proton_momentum_FluxUncertainty.toml',
    'flux_exp': f'{BASE}/leading_proton_momentum_FluxUncertainty_Exp.toml',
    'detsys':   f'{BASE}/leading_proton_momentum_DetectorUncertainty.toml',
}

def run(close_figs=True):
    for name, config in CONFIGS.items():
        print(f'\n--- Running: {name} ---')
        analysis = Analysis(config, DATA_FILE)
        analysis.run(close_figs=close_figs)

if __name__ == '__main__':
    run()
