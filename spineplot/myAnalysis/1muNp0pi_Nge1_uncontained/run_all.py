import os
import sys
import multiprocessing

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

print("All modules loaded. Starting parallel analysis...\n", flush=True)

OBSERVABLES = [
    ('dalphaT',              dalphaT_plots),
    ('dpT',                  dpT_plots),
    ('dphiT',                dphiT_plots),
    ('leading_muon_NuMI_angle',  leading_muon_NuMI_angle_plots),
    ('leading_muon_length',      leading_muon_length_plots),
    ('leading_muon_momentum',    leading_muon_momentum_plots),
    ('leading_muon_polar_angle', leading_muon_polar_angle_plots),
    ('leading_pion_length',      leading_pion_length_plots),
    ('leading_proton_momentum',  leading_proton_momentum_plots),
    ('opening_angle',            opening_angle_plots),
    ('vertex_x',                 vertex_x_plots),
    ('vertex_y',                 vertex_y_plots),
    ('vertex_z',                 vertex_z_plots),
]


def _run_one(args):
    name, module = args
    print(f"[{name}] Starting...", flush=True)
    module.run(close_figs=True)
    print(f"[{name}] Done.", flush=True)


def main():
    with multiprocessing.Pool(processes=len(OBSERVABLES)) as pool:
        pool.map(_run_one, OBSERVABLES)
    print("\nAll observables complete.", flush=True)


if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')
    main()
