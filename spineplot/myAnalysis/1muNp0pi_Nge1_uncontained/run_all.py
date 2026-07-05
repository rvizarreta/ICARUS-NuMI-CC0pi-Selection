import os
import sys
import multiprocessing
import importlib
import traceback

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
    'proton_multiplicity',
    'vertex_x',
    'vertex_y',
    'vertex_z',
]


def _run_one(name):
    try:
        sys.path.insert(0, os.path.join(workspace_root, '../../..'))
        sys.path.insert(0, os.path.join(workspace_root, name))
        print(f'[{name}] Starting...', flush=True)
        module = importlib.import_module(f'{name}_plots')
        module.run(close_figs=True)
        print(f'[{name}] Done.', flush=True)
    except Exception:
        print(f'[{name}] FAILED:\n{traceback.format_exc()}', flush=True)
        sys.exit(1)


N_PARALLEL = 4  # increase if GPVM has spare memory, reduce if it freezes


if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')
    print(f'Starting analysis: {len(OBSERVABLES)} observables, {N_PARALLEL} at a time...\n', flush=True)

    all_processes = [
        multiprocessing.Process(target=_run_one, args=(name,), name=name)
        for name in OBSERVABLES
    ]

    # Run in batches of N_PARALLEL
    for i in range(0, len(all_processes), N_PARALLEL):
        batch = all_processes[i:i + N_PARALLEL]
        names = [p.name for p in batch]
        print(f'--- Batch {i // N_PARALLEL + 1}: {names} ---', flush=True)
        for p in batch:
            p.start()
        for p in batch:
            p.join()

    failed = [p.name for p in all_processes if p.exitcode != 0]
    if failed:
        print(f'\nFailed observables: {failed}', flush=True)
    else:
        print('\nAll observables complete.', flush=True)
