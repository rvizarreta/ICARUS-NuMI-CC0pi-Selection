#%%
import toml
import os
import sys
workspace_root = os.getcwd()  
sys.path.insert(0, workspace_root + "/..")
from analysis import Analysis
#%%
# ANALYSIS SPECIFIC DIRECTORIES
ana_directory = "1muNp0pi_Nge1_uncontained"
#data_file = '/Users/rvizarreta/Library/CloudStorage/GoogleDrive-rvizarreta14@gmail.com/My Drive/🏛 PhD Repository/🚀 Research/🤖 Experiments&Projects/ICARUS/ICARUS_CC0pi_Selection/data/1muNp0pi_Nge1_uncontained.root'
data_file = '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root'
#data_file = '/Users/rvizarreta/Library/CloudStorage/GoogleDrive-rvizarreta14@gmail.com/My Drive/🏛 PhD Repository/🚀 Research/🤖 Experiments&Projects/ICARUS/ICARUS_CC0pi_Selection/data/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root'
#%%
# Visible energy plot
analysis = Analysis(ana_directory + '/visibleEnergy/visibleEnergy.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta p_T plots
analysis = Analysis(ana_directory + '/dpT/dpT.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta phiT plots
analysis = Analysis(ana_directory + '/dphiT/dphiT.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta alpha_T plots
analysis = Analysis(ana_directory + '/dalphaT/dalphaT.toml', data_file)
analysis.run(close_figs=False)
#%%
# Leading Muon Momentum
analysis = Analysis(ana_directory + '/leading_muon_momentum/leading_muon_momentum.toml', data_file)
analysis.run(close_figs=False)
#%%
# Leading Muon Length
analysis = Analysis(ana_directory + '/leading_muon_length/leading_muon_length.toml', data_file)
analysis.run(close_figs=False)
#%%
# Leading muon polar angle plots
analysis = Analysis(ana_directory + '/leading_muon_polar_angle/leading_muon_polar_angle.toml', data_file)
analysis.run(close_figs=False)
#%%
# Leading Proton Momentum
analysis = Analysis(ana_directory + '/leading_proton_momentum/leading_proton_momentum.toml', data_file)
analysis.run(close_figs=False)
#%%
# Opening Angle
analysis = Analysis(ana_directory + '/opening_angle/opening_angle.toml', data_file)
analysis.run(close_figs=False)
#%%
# Vertex X
analysis = Analysis(ana_directory + '/vertex_x/vertex_x.toml', data_file)
analysis.run(close_figs=False)
#%%
# Vertex Y
analysis = Analysis(ana_directory + '/vertex_y/vertex_y.toml', data_file)
analysis.run(close_figs=False)
#%%
# Vertex Z
analysis = Analysis(ana_directory + '/vertex_z/vertex_z.toml', data_file)
analysis.run(close_figs=False)
#%%
analysis = Analysis(ana_directory + '/efficiency.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta p_T cross-section uncertainties plot
analysis = Analysis(ana_directory + '/dpT/dpT_xSecUncertainty.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta p_T cross-section uncertainties plot
analysis = Analysis(ana_directory + '/dpT/dpT_xSecUncertainty_Exp.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta p_T cross-section uncertainties plot
analysis = Analysis(ana_directory + '/dpT/dpT_FluxUncertainty.toml', data_file)
analysis.run(close_figs=False)
#%%
# Delta p_T cross-section uncertainties plot
analysis = Analysis(ana_directory + '/dpT/dpT_FluxUncertainty_Exp.toml', data_file)
analysis.run(close_figs=False)
#%%
# No charged pions - Truth
analysis = Analysis(ana_directory + '/pions_in_otherCC/pions_in_otherCC.toml', data_file)
analysis.run(close_figs=False)
#%%
