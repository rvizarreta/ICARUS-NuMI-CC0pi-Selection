#!/bin/bash
LOG=/tmp/plotjob.$$.log
exec > "$LOG" 2>&1
set -x
PROJECT="$1"

OBSERVABLES=(
    dalphaT dpT dphiT
    leading_muon_NuMI_angle leading_muon_length leading_muon_momentum
    leading_muon_polar_angle leading_pion_length leading_proton_momentum
    opening_angle vertex_x vertex_y vertex_z
)
NANALYSES=6

OBS_INDEX=$(( PROCESS / NANALYSES ))
ANALYSIS_INDEX=$(( PROCESS % NANALYSES ))
OBS=${OBSERVABLES[$OBS_INDEX]}
TAG="${OBS}_a${ANALYSIS_INDEX}"

source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup sbnana v10_01_02_01 -q e26:prof
setup cmake v3_27_4
PY=$(which python3)
export MPLBACKEND=Agg
export MPLCONFIGDIR=$_CONDOR_SCRATCH_DIR/mpl
mkdir -p "$MPLCONFIGDIR"

STATUS=$PROJECT/status
mark() {
  local f=/tmp/mark.$$
  echo "$(date +%H:%M:%S) PROC=$PROCESS OBS=$OBS AIDX=$ANALYSIS_INDEX step=$1" > "$f"
  ifdh cp --cp_maxretries=0 --web_timeout=100 "$f" "$STATUS/${TAG}_$1"
}
ship_log() {
  ifdh cp --cp_maxretries=0 --web_timeout=100 "$LOG" "$PROJECT/output/run_${CLUSTER}_${TAG}.log"
}
trap ship_log EXIT

mark 00_start
[ -z "$OBS" ] && { echo "OBS empty for PROCESS=$PROCESS"; mark 99_bad_obs; exit 1; }

cd "$_CONDOR_SCRATCH_DIR" || exit 1

export PYTHONUSERBASE=$_CONDOR_SCRATCH_DIR/pylocal
mkdir -p "$PYTHONUSERBASE"
mark 01_pip_start
"$PY" -m pip install --user --quiet \
    uproot awkward numpy matplotlib toml \
    hist boost-histogram mplhep \
    pandas scipy tqdm scikit-learn
PIP_RC=$?
echo "pip rc=$PIP_RC"
[ $PIP_RC -ne 0 ] && { mark 99_pip_failed; exit 1; }
mark 02_pip_done

git clone https://github.com/rvizarreta/ICARUS-NuMI-CC0pi-Selection.git
cd ICARUS-NuMI-CC0pi-Selection || { mark 99_clone_failed; exit 1; }
git checkout feature/rvizarr_cc0pi_selection
echo "checkout rc=$?"
mark 03_clone_done

REPO_DIR=$(pwd)
ANA_DIR=$REPO_DIR/spineplot/myAnalysis/1muNp0pi_Nge1_uncontained

LOCAL_DATA=$_CONDOR_SCRATCH_DIR/input.root
mark 04_stage_start
ifdh cp --cp_maxretries=0 --web_timeout=100 \
    /pnfs/icarus/persistent/users/rvizarr/plots/icarus_numi_numu_mc_onbeam_offbeam_syst_ppfx.root \
    "$LOCAL_DATA"
STAGE_RC=$?
echo "stage rc=$STAGE_RC"
[ $STAGE_RC -ne 0 ] && { mark 99_stage_failed; exit 1; }
[ -s "$LOCAL_DATA" ] || { echo "staged file empty/missing"; mark 99_stage_empty; exit 1; }
export SPINE_DATA_FILE=$LOCAL_DATA
mark 05_stage_done

LOCAL_G4=$_CONDOR_SCRATCH_DIR/g4
mkdir -p "$LOCAL_G4"

mark 05a_g4_stage_start
for f in /pnfs/icarus/persistent/users/rvizarr/plots/G4/*.root; do
    ifdh cp --cp_maxretries=0 --web_timeout=100 "$f" "$LOCAL_G4/$(basename "$f")"
done
G4_COUNT=$(ls -1 "$LOCAL_G4" 2>/dev/null | wc -l)
echo "staged $G4_COUNT G4 root files"
[ "$G4_COUNT" -eq 30 ] || { echo "expected 30 G4 root files, got $G4_COUNT"; mark 99_g4_stage_failed; exit 1; }

find "$ANA_DIR" -name '*.toml' -exec sed -i \
    "s|/exp/icarus/app/users/rvizarr/gundam-icarus/configs/Configs_ParameterSet/G4/outputs|$LOCAL_G4|g" \
    {} +
mark 05b_g4_stage_done

cd "$ANA_DIR/$OBS" || { mark 99_cd_failed; exit 1; }
touch /tmp/before_plot.$$
mark 06_plot_start
"$PY" -c "
import sys
sys.path.insert(0, '$REPO_DIR/spineplot')
sys.path.insert(0, '$ANA_DIR/$OBS')
import importlib
m = importlib.import_module('${OBS}_plots')
m.run_one(m.ANALYSES[$ANALYSIS_INDEX])
"
PLOT_RC=$?
echo "python rc=$PLOT_RC"
if [ $PLOT_RC -ne 0 ]; then
    mark 99_plot_failed
    exit 1
fi
mark 07_plot_done

OUT_JPEG="$ANA_DIR/$OBS/jpeg"
OUT_PDF="$ANA_DIR/$OBS/pdf"
echo "=== jpeg dir ==="; ls -la "$OUT_JPEG" 2>&1
echo "=== pdf dir ===";  ls -la "$OUT_PDF"  2>&1

shipped=0
shopt -s nullglob
for f in "$OUT_JPEG"/*.jpeg "$OUT_PDF"/*.pdf; do
    [ "$f" -nt /tmp/before_plot.$$ ] || continue
    ifdh cp --cp_maxretries=0 --web_timeout=100 "$f" "$PROJECT/output/${TAG}__$(basename "$f")"
    echo "ship $(basename "$f") rc=$?"
    shipped=$((shipped+1))
done
echo "shipped $shipped files"
if [ $shipped -eq 0 ]; then
    mark 99_no_output
    exit 1
fi
mark 08_ship_done
