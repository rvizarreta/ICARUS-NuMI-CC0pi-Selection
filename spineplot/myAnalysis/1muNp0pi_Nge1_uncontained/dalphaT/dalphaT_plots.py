import toml
import os
import sys

workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))
from analysis import Analysis

data_file = '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root'

ana_directory = os.path.dirname(os.path.abspath(__file__))

def run(close_figs=True):
    Analysis(os.path.join(ana_directory, 'dalphaT.toml'), data_file).run(close_figs=close_figs)
    Analysis(os.path.join(ana_directory, 'dalphaT_xSecUncertainty.toml'), data_file).run(close_figs=close_figs)
    Analysis(os.path.join(ana_directory, 'dalphaT_xSecUncertainty_Exp.toml'), data_file).run(close_figs=close_figs)
    Analysis(os.path.join(ana_directory, 'dalphaT_FluxUncertainty.toml'), data_file).run(close_figs=close_figs)
    Analysis(os.path.join(ana_directory, 'dalphaT_FluxUncertainty_Exp.toml'), data_file).run(close_figs=close_figs)
    Analysis(os.path.join(ana_directory, 'dalphaT_DetectorUncertainty.toml'), data_file).run(close_figs=close_figs)

if __name__ == '__main__':
    run()