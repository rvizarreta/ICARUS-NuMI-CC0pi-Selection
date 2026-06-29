import os
import sys

workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))

sys.path.insert(0, os.path.join(workspace_root, "dalphaT"))
import dalphaT_plots

sys.path.insert(0, os.path.join(workspace_root, "dpT"))
import dpT_plots

sys.path.insert(0, os.path.join(workspace_root, "dphiT"))
import dphiT_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_NuMI_angle"))
import leading_muon_NuMI_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_length"))
import leading_muon_length_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_momentum"))
import leading_muon_momentum_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_polar_angle"))
import leading_muon_polar_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_pion_length"))
import leading_pion_length_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_proton_momentum"))
import leading_proton_momentum_plots

sys.path.insert(0, os.path.join(workspace_root, "opening_angle"))
import opening_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_x"))
import vertex_x_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_y"))
import vertex_y_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_z"))
import vertex_z_plots


def main(close_figs=True):
    print("=== Running dalphaT plots ===")
    dalphaT_plots.run(close_figs=close_figs)

    print("=== Running dpT plots ===")
    dpT_plots.run(close_figs=close_figs)

    print("=== Running dphiT plots ===")
    dphiT_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_NuMI_angle plots ===")
    leading_muon_NuMI_angle_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_length plots ===")
    leading_muon_length_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_momentum plots ===")
    leading_muon_momentum_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_polar_angle plots ===")
    leading_muon_polar_angle_plots.run(close_figs=close_figs)

    print("=== Running leading_pion_length plots ===")
    leading_pion_length_plots.run(close_figs=close_figs)

    print("=== Running leading_proton_momentum plots ===")
    leading_proton_momentum_plots.run(close_figs=close_figs)

    print("=== Running opening_angle plots ===")
    opening_angle_plots.run(close_figs=close_figs)

    print("=== Running vertex_x plots ===")
    vertex_x_plots.run(close_figs=close_figs)

    print("=== Running vertex_y plots ===")
    vertex_y_plots.run(close_figs=close_figs)

    print("=== Running vertex_z plots ===")
    vertex_z_plots.run(close_figs=close_figs)


if __name__ == '__main__':
    main()
