import os
import sys

workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))
from analysis import Analysis

data_file = os.environ.get('SPINE_DATA_FILE', '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root')

ana_directory = os.path.dirname(os.path.abspath(__file__))

ANALYSES = [
    'leading_muon_polar_angle.toml',
    'leading_muon_polar_angle_xSecUncertainty.toml',
    'leading_muon_polar_angle_xSecUncertainty_Exp.toml',
    'leading_muon_polar_angle_FluxUncertainty.toml',
    'leading_muon_polar_angle_FluxUncertainty_Exp.toml',
    'leading_muon_polar_angle_DetectorUncertainty.toml',
]

def run_one(toml_name, close_figs=True):
    Analysis(os.path.join(ana_directory, toml_name), data_file).run(close_figs=close_figs)

def run(close_figs=True):
    for toml_name in ANALYSES:
        run_one(toml_name, close_figs=close_figs)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_one(ANALYSES[int(sys.argv[1])])
    else:
        run()
