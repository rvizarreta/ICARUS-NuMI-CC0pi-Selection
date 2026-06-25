#!/usr/bin/env python
"""
generate_NormCCRES_dial.py

Add a synthetic NormCCRES dial branch to the GUNDAM input ROOT file.

For each event, stores a TClonesArray containing one TGraph, matching the
format of the existing GENIEReWeight_SBN_v1_multisigma_* branches.

  For RES events (mode == 1): TGraph values are (1 + σ × SIGMA_WIDTH),
                              σ ∈ [-3, -2, -1, 0, +1, +2, +3].
  For all other events:        TGraph values are 1.0 at every σ.

The CC νμ selection is already enforced upstream (signal/sideband selections),
so no explicit CC cut is applied here.
"""
import ROOT
import array

INPUT_FILE  = (
    "/Users/rvizarreta/Library/CloudStorage/"
    "GoogleDrive-rvizarreta14@gmail.com/My Drive/"
    "🏛 PhD Repository/🚀 Research/🤖 Experiments&Projects/"
    "ICARUS/ICARUS_CC0pi_Selection/data/"
    "icarus_numi_numu_mc_onbeam_offbeam_syst_gundam_withPCA.root"
)
OUTPUT_FILE = (
    "/Users/rvizarreta/Library/CloudStorage/"
    "GoogleDrive-rvizarreta14@gmail.com/My Drive/"
    "🏛 PhD Repository/🚀 Research/🤖 Experiments&Projects/"
    "ICARUS/ICARUS_CC0pi_Selection/data/"
    "icarus_numi_numu_mc_onbeam_offbeam_syst_gundam_withPCA_NormCCRES.root"
)

SIGMA_WIDTH = 0.20   # 20% per σ — placeholder, refine from generator spread

SIGMA_POINTS   = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
N_SIGMA        = len(SIGMA_POINTS)
RES_MODE_VALUE = 1   # NUISANCE convention in this file: mode == 1 ⇒ RES

MC_TREES = [
    "events/full/selected",
    "events/full/sideband",
    "events/full/signal",
]
COPY_OBJECTS = [
    "events/full/POT",
    "events/full/Livetime",
    "events/onbeam/selected",
    "events/onbeam/sideband",
    "events/onbeam/POT",
    "events/onbeam/Livetime",
]


def get_or_make_dir(rfile, parts):
    d = rfile
    for p in parts:
        sub = d.GetDirectory(p)
        if not sub:
            sub = d.mkdir(p)
        d = sub
    return d


def copy_object_as_is(f_in, f_out, path):
    obj = f_in.Get(path)
    if not obj:
        print(f"  [skip] {path} not found")
        return
    parts = path.split("/")
    d = get_or_make_dir(f_out, parts[:-1])
    d.cd()

    if isinstance(obj, ROOT.TTree):
        out = obj.CloneTree(-1, "fast")
        out.Write()
        print(f"  [copy] {path}: TTree, {out.GetEntries()} entries")
    else:
        # Histograms, graphs, anything non-tree
        out = obj.Clone()
        if hasattr(out, "SetDirectory"):
            out.SetDirectory(d)
        out.Write()
        print(f"  [copy] {path}: {type(out).__name__}")


def add_dial_branch(f_in, f_out, path):
    t_in = f_in.Get(path)
    if not t_in:
        print(f"  [skip] {path} not found")
        return
    parts = path.split("/")
    get_or_make_dir(f_out, parts[:-1]).cd()

    t_out = t_in.CloneTree(0)

    tca = ROOT.TClonesArray("TGraph", 1)
    t_out.Branch("Synthetic_multisigma_NormCCRES",
                 "TClonesArray", ROOT.AddressOf(tca), 32000, 0)

    fX = array.array('d', SIGMA_POINTS)
    fY_active   = array.array('d', [1.0 + s * SIGMA_WIDTH for s in SIGMA_POINTS])
    fY_inactive = array.array('d', [1.0] * N_SIGMA)

    n_total = t_in.GetEntries()
    n_res   = 0

    for i in range(n_total):
        t_in.GetEntry(i)
        is_res = (t_in.true_interaction_mode == RES_MODE_VALUE)
        if is_res:
            n_res += 1
            fY = fY_active
        else:
            fY = fY_inactive

        tca.Clear()
        tca.ExpandCreate(1)
        tg = tca.At(0)
        tg.Set(N_SIGMA)
        for j in range(N_SIGMA):
            tg.SetPoint(j, fX[j], fY[j])

        t_out.Fill()

        if i % 20000 == 0 and i > 0:
            print(f"    {i}/{n_total}", flush=True)

    t_out.Write()
    pct = 100.0 * n_res / max(n_total, 1)
    print(f"  [add ] {path}: {n_total} entries, "
          f"{n_res} RES ({pct:.1f}%)")


def main():
    print(f"Input:  {INPUT_FILE}")
    print(f"Output: {OUTPUT_FILE}")
    print(f"σ width: ±{100*SIGMA_WIDTH:.0f}% per σ")

    f_in  = ROOT.TFile.Open(INPUT_FILE, "READ")
    f_out = ROOT.TFile(OUTPUT_FILE, "RECREATE")

    print("\nCopying ancillary objects:")
    for path in COPY_OBJECTS:
        copy_object_as_is(f_in, f_out, path)

    print("\nProcessing MC trees:")
    for path in MC_TREES:
        add_dial_branch(f_in, f_out, path)

    f_out.Close()
    f_in.Close()
    print("\nDone.")


if __name__ == "__main__":
    main()