import os
import sys

workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
sys.path.insert(0, workspace_root)

from analysis import Analysis

DATA_FILE = '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root'

BASE = os.path.dirname(os.path.abspath(__file__))

ANALYSES = [
    'leading_proton_momentum.toml',
    'leading_proton_momentum_xSecUncertainty.toml',
    'leading_proton_momentum_xSecUncertainty_Exp.toml',
    'leading_proton_momentum_FluxUncertainty.toml',
    'leading_proton_momentum_FluxUncertainty_Exp.toml',
    'leading_proton_momentum_DetectorUncertainty.toml',
]

def run_one(toml_name, close_figs=True):
    Analysis(os.path.join(BASE, toml_name), DATA_FILE).run(close_figs=close_figs)

def run(close_figs=True):
    for toml_name in ANALYSES:
        run_one(toml_name, close_figs=close_figs)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_one(ANALYSES[int(sys.argv[1])])
    else:
        run()
