import os
import sys
import multiprocessing
import importlib

workspace_root = os.path.dirname(os.path.abspath(__file__))

OBSERVABLES = [
    'dalphaT',
    'dpT',
    'dphiT',
    'leading_muon_NuMI_angle',
    'leading_muon_length',
    'leading_muon_momentum',
    'leading_muon_polar_angle',
    'leading_pion_length',
    'leading_proton_momentum',
    'opening_angle',
    'vertex_x',
    'vertex_y',
    'vertex_z',
]


def _run_one(name):
    sys.path.insert(0, os.path.join(workspace_root, '../../..'))
    sys.path.insert(0, os.path.join(workspace_root, name))
    print(f'[{name}] Starting...', flush=True)
    module = importlib.import_module(f'{name}_plots')
    module.run(close_figs=True)
    print(f'[{name}] Done.', flush=True)


if __name__ == '__main__':
    print(f'Starting parallel analysis for {len(OBSERVABLES)} observables...\n', flush=True)
    multiprocessing.set_start_method('spawn')
    with multiprocessing.Pool(processes=len(OBSERVABLES)) as pool:
        pool.map(_run_one, OBSERVABLES)
    print('\nAll observables complete.', flush=True)
