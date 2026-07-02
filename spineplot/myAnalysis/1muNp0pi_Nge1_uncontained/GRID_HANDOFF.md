# Grid Submission Handoff

## What has been done

All 13 `<observable>_plots.py` files in `spineplot/myAnalysis/1muNp0pi_Nge1_uncontained/` have been refactored and standardized. Each file now has:

- A uniform header (same imports, `workspace_root`, `sys.path`, `ana_directory`)
- `data_file = os.environ.get('SPINE_DATA_FILE', '/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root')`  
  → reads from env var on the grid, falls back to `/exp/` path for local/GPVM use
- `ANALYSES` list: 6 TOML filenames per observable (main + xSec + xSec_Exp + Flux + Flux_Exp + DetSys)
- `run_one(toml_name, close_figs=True)`: runs a single TOML config
- `run(close_figs=True)`: runs all 6 sequentially (unchanged behavior)
- `__main__` block: `python3 <obs>_plots.py 0` runs ANALYSES[0], no arg runs all

## Observables and their ANALYSES indices

| Index | TOML suffix |
|-------|-------------|
| 0 | `<obs>.toml` (main) |
| 1 | `<obs>_xSecUncertainty.toml` |
| 2 | `<obs>_xSecUncertainty_Exp.toml` |
| 3 | `<obs>_FluxUncertainty.toml` |
| 4 | `<obs>_FluxUncertainty_Exp.toml` |
| 5 | `<obs>_DetectorUncertainty.toml` |

Observables (13 total):
`dalphaT`, `dpT`, `dphiT`, `leading_muon_NuMI_angle`, `leading_muon_length`, `leading_muon_momentum`, `leading_muon_polar_angle`, `leading_pion_length`, `leading_proton_momentum`, `opening_angle`, `vertex_x`, `vertex_y`, `vertex_z`

Total grid jobs: **13 × 6 = 78**

## Key paths (GPVM)

- **Repo**: `/exp/icarus/app/users/rvizarr/ICARUS-NuMI-CC0pi-Selection/`
- **Spineplot dir**: `<repo>/spineplot/`
- **Analysis dir**: `<repo>/spineplot/myAnalysis/1muNp0pi_Nge1_uncontained/`
- **Data file (current)**: `/exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root`
- **Conda env**: `plotting_env` located at `<repo>/spineplot/myAnalysis/plotting_env/`
- **Output**: plots saved to `<observable>/jpeg/` and `<observable>/pdf/` relative to each TOML's directory

## Next step: grid submission

Worker nodes cannot access `/exp/`. To run on the grid we need to:

1. **Stage data to `/pnfs/`**  
   ```bash
   ifdh cp /exp/icarus/data/users/rvizarr/medulla/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root \
       /pnfs/icarus/scratch/users/rvizarr/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root
   ```

2. **Make the conda env available on worker nodes**  
   Options: tarball it and stage to `/pnfs/`, or install a fresh env inside the job.

3. **Write a wrapper shell script** that:
   - Sets up the environment
   - Sets `SPINE_DATA_FILE` to the `/pnfs/` path
   - Calls `python3 <obs>/<obs>_plots.py <index>`
   - Stages output plots to `/pnfs/`

4. **Submit 78 jobs** (one per observable × TOML index) via `jobsub_submit`

## Git branch

`feature/rvizarr_cc0pi_selection` — all changes committed and pushed.
