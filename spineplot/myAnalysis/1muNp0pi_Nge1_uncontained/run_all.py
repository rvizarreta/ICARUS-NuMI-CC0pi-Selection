import os
import sys

workspace_root = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(workspace_root, "../../.."))

print("Importing modules...", flush=True)

sys.path.insert(0, os.path.join(workspace_root, "dalphaT"))
print("  dalphaT", flush=True)
import dalphaT_plots

sys.path.insert(0, os.path.join(workspace_root, "dpT"))
print("  dpT", flush=True)
import dpT_plots

sys.path.insert(0, os.path.join(workspace_root, "dphiT"))
print("  dphiT", flush=True)
import dphiT_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_NuMI_angle"))
print("  leading_muon_NuMI_angle", flush=True)
import leading_muon_NuMI_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_length"))
print("  leading_muon_length", flush=True)
import leading_muon_length_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_momentum"))
print("  leading_muon_momentum", flush=True)
import leading_muon_momentum_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_muon_polar_angle"))
print("  leading_muon_polar_angle", flush=True)
import leading_muon_polar_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_pion_length"))
print("  leading_pion_length", flush=True)
import leading_pion_length_plots

sys.path.insert(0, os.path.join(workspace_root, "leading_proton_momentum"))
print("  leading_proton_momentum", flush=True)
import leading_proton_momentum_plots

sys.path.insert(0, os.path.join(workspace_root, "opening_angle"))
print("  opening_angle", flush=True)
import opening_angle_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_x"))
print("  vertex_x", flush=True)
import vertex_x_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_y"))
print("  vertex_y", flush=True)
import vertex_y_plots

sys.path.insert(0, os.path.join(workspace_root, "vertex_z"))
print("  vertex_z", flush=True)
import vertex_z_plots

print("All modules loaded. Starting analysis...\n", flush=True)


def main(close_figs=True):
    print("=== Running dalphaT plots ===", flush=True)
    dalphaT_plots.run(close_figs=close_figs)

    print("=== Running dpT plots ===", flush=True)
    dpT_plots.run(close_figs=close_figs)

    print("=== Running dphiT plots ===", flush=True)
    dphiT_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_NuMI_angle plots ===", flush=True)
    leading_muon_NuMI_angle_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_length plots ===", flush=True)
    leading_muon_length_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_momentum plots ===", flush=True)
    leading_muon_momentum_plots.run(close_figs=close_figs)

    print("=== Running leading_muon_polar_angle plots ===", flush=True)
    leading_muon_polar_angle_plots.run(close_figs=close_figs)

    print("=== Running leading_pion_length plots ===", flush=True)
    leading_pion_length_plots.run(close_figs=close_figs)

    print("=== Running leading_proton_momentum plots ===", flush=True)
    leading_proton_momentum_plots.run(close_figs=close_figs)

    print("=== Running opening_angle plots ===", flush=True)
    opening_angle_plots.run(close_figs=close_figs)

    print("=== Running vertex_x plots ===", flush=True)
    vertex_x_plots.run(close_figs=close_figs)

    print("=== Running vertex_y plots ===", flush=True)
    vertex_y_plots.run(close_figs=close_figs)

    print("=== Running vertex_z plots ===", flush=True)
    vertex_z_plots.run(close_figs=close_figs)


if __name__ == '__main__':
    main()
