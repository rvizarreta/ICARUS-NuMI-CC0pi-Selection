#!/usr/bin/env python3
"""
NuMI spill-quality cuts (python port of B. Howard's SpillQualityCuts.cxx / BeamExposureVars.cxx,
stringent mode) applied to the SPINE on-beam flat CAFs.

For every file it produces
  * an EVENT table  (run, subrun, evt, cut flags) using hdr.spillnumiinfo (the triggering spill)
  * a RECORD table  (run, subrun, event, pot, cut flags) using hdr.numiinfo (all spills of the subrun,
                     stored only on the first event of each subrun)
and prints the POT after each cut.

Cut flag bits: 1 = POT in spill, 2 = horn current, 4 = beam position at target, 8 = beam width
(the width fit is only run on spills that already pass the first three, bit 8 is 0 otherwise).
A spill passes everything if flags == 15.
"""
import argparse, glob, os, sys, time
import numpy as np
import awkward as ak
import uproot
from multiprocessing import Pool
from scipy.optimize import curve_fit

P = "rec.hdr."
SCALARS = ["NSLINA", "NSLINB", "NSLINC", "NSLIND", "TRTGTD", "TR101D"]
ARR7 = ["HP121", "VP121", "HPTGT", "HITGT", "VPTGT", "VITGT"]
NMTGT = 199  # MTGTDS needs at least 199 channels (Howard: size < 199 -> fail)

MM_FT = 12. * 2.54 * 10.
Z_HP121, Z_VP121 = -68.04458 * MM_FT, -66.99283 * MM_FT
Z_HPTGT, Z_VPTGT = -31.25508 * MM_FT, -30.16533 * MM_FT
XPOS = 0.5 * (np.arange(48) + 1 - 24.5)


def extrap(v1, l1, v2, l2, l3):
    return v1 + (l3 - l1) * (v2 - v1) / (l2 - l1)


# ---------------------------------------------------------------- cuts (vectorised)
def pot_val(tg, tr):
    return np.where(tg > 0.02, tg, np.where(tr >= 0., tr, 0.))


def horn_val(a, b, c, d):
    return (a - 0.01) / 0.9951 + (b + 0.14) / 0.9957 + (c + 0.05) / 0.9965 + (d + 0.07) / 0.9945


def beam_pos(run, hp121, vp121, hptgt, hitgt, vptgt, vitgt, ok):
    """beam position at the target (x,y) in mm; NaN where the monitors are unusable"""
    wh, wv = hitgt.sum(1), vitgt.sum(1)
    good = ok & (wh > np.finfo(float).eps) & (wv > np.finfo(float).eps)
    with np.errstate(all="ignore"):
        ave_hp121, ave_vp121 = hp121.mean(1), vp121.mean(1)
        ave_hptgt = (hptgt * hitgt).sum(1) / wh
        ave_vptgt = (vptgt * vitgt).sum(1) / wv
    # offsets (Run1 and earlier: run <= 8553, Run2: run > 8553). NB: Howard's yCorr uses the H z-positions too.
    r1 = run <= 8553
    x121 = np.where(r1, 0., 1.2); xtgt = np.where(r1, 0., 0.7)
    y121 = np.where(r1, -1.5, -0.4); ytgt = np.where(r1, -1.5, -0.94)
    xnom = np.where(r1, 0.398, 0.03); ynom = np.where(r1, -0.39, -0.59)
    xcorr = extrap(x121, Z_HP121, xtgt, Z_HPTGT, 0.) - xnom
    ycorr = extrap(y121, Z_HP121, ytgt, Z_HPTGT, 0.) - ynom
    x = extrap(ave_hp121, Z_HP121, ave_hptgt, Z_HPTGT, 0.) - xcorr
    y = extrap(ave_vp121, Z_VP121, ave_vptgt, Z_VPTGT, 0.) - ycorr
    x = np.where(good, x, np.nan); y = np.where(good, y, np.nan)
    return x, y


def _gauss(x, p0, p1, p2, p3):
    return p0 + p1 * np.exp(-((x - p2) * (x - p2)) / (2. * p3 * p3))


def sigma_fit(chan, mode="howard"):
    """sigma of Gaussian+offset fit, same model / sigma limit as Howard's TF1 fit.
    mode 'howard': Howard's start values (p0=-mean, p1=-min, mu=0, sigma=2)
    mode 'data'  : start values taken from the profile itself (baseline = median, amplitude = largest excursion)"""
    if mode == "howard":
        p0 = [-chan.sum() / 48., -chan.min(), 0., 2.]
    else:
        base = np.median(chan); j = np.argmax(np.abs(chan - base))
        p0 = [base, chan[j] - base, XPOS[j], 2.]
    try:
        popt, _ = curve_fit(_gauss, XPOS, chan, p0=p0,
                            bounds=([-np.inf, -np.inf, -np.inf, 0.], [np.inf, np.inf, np.inf, 10.]),
                            maxfev=2000)
        return popt[3]
    except Exception:
        return np.nan


def beam_width(mtgtds, ok, mode="howard"):
    n = len(mtgtds)
    sh = np.full(n, np.nan); sv = np.full(n, np.nan)
    for i in np.nonzero(ok)[0]:
        sh[i] = sigma_fit(-mtgtds[i, 103:151], mode)
        sv[i] = sigma_fit(-mtgtds[i, 151:199], mode)
    return sh, sv


def apply_cuts(run, d, do_width=True, mode="howard"):
    """d: dict with scalars (N,), arrays (N,7) / mtgtds (N,199) and *_ok masks. Returns pot, flags."""
    pot = pot_val(d["TRTGTD"], d["TR101D"])
    f_pot = pot > 2.0e12
    horn = horn_val(d["NSLINA"], d["NSLINB"], d["NSLINC"], d["NSLIND"])
    f_horn = (horn > -202.) & (horn <= -196.4)
    x, y = beam_pos(run, d["HP121"], d["VP121"], d["HPTGT"], d["HITGT"], d["VPTGT"], d["VITGT"], d["ok7"])
    with np.errstate(invalid="ignore"):
        f_pos = (np.abs(x) <= 1.) & (np.abs(y) <= 1.)
    flags = f_pot.astype(np.uint8) | (f_horn.astype(np.uint8) << 1) | (f_pos.astype(np.uint8) << 2)
    if do_width:
        todo = (flags == 7) & d["okM"]
        sh, sv = beam_width(d["MTGTDS"], todo, mode)
        with np.errstate(invalid="ignore"):
            f_w = (sh > 0.57) & (sh <= 1.88) & (sv > 0.57) & (sv <= 1.88)
        flags = flags | (f_w.astype(np.uint8) << 3)
    return pot, flags


# ---------------------------------------------------------------- reading
def to2d(flat, lens, L):
    """flat values + per-row lengths -> (rows, L) array (NaN padded / clipped) and mask lens>=L"""
    a = ak.unflatten(flat, lens)
    a = ak.fill_none(ak.pad_none(a, L, clip=True), np.nan)
    return ak.to_numpy(a), lens >= L


def build(prefix, arr, per_record, nnumi=None):
    """collect the fields of hdr.<prefix> into plain numpy dicts"""
    def get(name):
        x = arr[P + prefix + "." + name]
        return ak.to_numpy(ak.flatten(x, axis=None)) if per_record else ak.to_numpy(x)
    d = {k: get(k).astype(np.float64) for k in SCALARS}
    ok7 = None
    for k in ARR7 + ["MTGTDS"]:
        x = arr[P + prefix + "." + k]; l = arr[P + prefix + "." + k + "..length"]
        flat = ak.to_numpy(ak.flatten(x, axis=None)).astype(np.float64)
        lens = ak.to_numpy(ak.flatten(l, axis=None) if per_record else l).astype(np.int64)
        L = NMTGT if k == "MTGTDS" else 7
        d[k], m = to2d(flat, lens, L)
        if k == "MTGTDS":
            d["okM"] = m
        else:
            m7 = lens == 7
            ok7 = m7 if ok7 is None else (ok7 & m7)
    d["ok7"] = ok7
    d["event"] = get("event").astype(np.int64)
    d["time"] = get("spill_time_s").astype(np.float64) + 1e-9 * get("spill_time_ns").astype(np.float64)
    return d


def process(args):
    path, do_width, mode = args
    t0 = time.time()
    try:
        f = uproot.open(path)
        tn = [k for k in f.keys() if "recTree" in k or "Events" in k][0]
        t = f[tn]
        names = [P + "run", P + "subrun", P + "evt", P + "pot", P + "nnumiinfo"]
        for pre in ("numiinfo", "spillnumiinfo"):
            names += [P + pre + "." + k for k in SCALARS + ["event", "spill_time_s", "spill_time_ns"]]
            names += [P + pre + "." + k for k in ARR7 + ["MTGTDS"]]
            names += [P + pre + "." + k + "..length" for k in ARR7 + ["MTGTDS"]]
        a = t.arrays(names, library="ak")
        run = ak.to_numpy(a[P + "run"]).astype(np.int64)
        sub = ak.to_numpy(a[P + "subrun"]).astype(np.int64)
        evt = ak.to_numpy(a[P + "evt"]).astype(np.int64)
        pot_hdr = float(ak.to_numpy(a[P + "pot"]).astype(np.float64).sum())
        nn = ak.to_numpy(a[P + "nnumiinfo"]).astype(np.int64)

        # --- triggering spill, one per entry
        ds = build("spillnumiinfo", a, False)
        spot, sflags = apply_cuts(run, ds, do_width, mode)
        ev = dict(run=run, subrun=sub, evt=evt, spill_event=ds["event"], pot=spot, flags=sflags, time=ds["time"])

        # --- all spill records of the subrun (first event only)
        rrun = np.repeat(run, nn); rsub = np.repeat(sub, nn)
        dr = build("numiinfo", a, True)
        assert len(dr["event"]) == nn.sum(), "numiinfo length mismatch in " + path
        rpot, rflags = apply_cuts(rrun, dr, do_width, mode)
        rec = dict(run=rrun, subrun=rsub, event=dr["event"], pot=rpot, flags=rflags, time=dr["time"])
        return path, ev, rec, pot_hdr, time.time() - t0, None
    except Exception as e:
        return path, None, None, 0., time.time() - t0, repr(e)


def key(run, sub, evt):
    return (run.astype(np.int64) << 40) | (sub.astype(np.int64) << 26) | evt.astype(np.int64)


def run_block(files, o):
    jobs = [(f, not o.no_width, o.fit_init) for f in files]
    res = Pool(o.nproc).map(process, jobs) if o.nproc > 1 else [process(j) for j in jobs]
    bad = [r[0] for r in res if r[5]]
    for r in res:
        if r[5]: print("FAILED", os.path.basename(r[0]), r[5], flush=True)
    ok = [r for r in res if not r[5]]
    if not ok:
        return None
    cat = lambda dicts, k: np.concatenate([d[k] for d in dicts])
    ev = {k: cat([r[1] for r in ok], k) for k in ok[0][1]}
    rec = {k: cat([r[2] for r in ok], k) for k in ok[0][2]}
    return ev, rec, sum(r[3] for r in ok), np.array([r[4] for r in ok]), bad


def save_npz(path, ev, rec, pot_hdr, secs, bad):
    np.savez(path, **{"ev_" + k: v for k, v in ev.items()}, **{"rec_" + k: v for k, v in rec.items()},
             pot_hdr=pot_hdr, secs=secs, bad=np.array(bad, dtype=str))


def load_npz(path):
    d = np.load(path)
    ev = {k[3:]: d[k] for k in d.files if k.startswith("ev_")}
    rec = {k[4:]: d[k] for k in d.files if k.startswith("rec_")}
    return ev, rec, float(d["pot_hdr"]), d["secs"], list(d["bad"])


def summarize(ev, rec, pot_hdr, secs, bad, good_runs=None, nfiles=None):
    kev = key(ev["run"], ev["subrun"], ev["evt"])
    print("\nfiles read OK: %d | failed: %d" % (nfiles if nfiles is not None else -1, len(bad)))
    for b in bad: print("   FAILED:", b)
    print("duplicate event keys:", len(kev) - len(np.unique(kev)))
    print("mean seconds per file: %.2f" % secs.mean())
    print("events: %d | triggering spill passes all cuts: %d (%.2f%%)" %
          (len(kev), (ev["flags"] == 15).sum(), 100. * (ev["flags"] == 15).mean()))
    for name, bit in (("POT", 1), ("horn", 2), ("position", 4), ("width", 8)):
        print("   event fails %-8s: %d" % (name, (((ev["flags"] & bit) == 0)).sum()))
    krec = key(rec["run"], rec["subrun"], rec["event"])
    inev = np.isin(krec, kev)
    print("\nspill records: %d | with event present in these files: %d" % (len(krec), inev.sum()))
    print("POT sum of hdr.pot (uncut, first_in_subrun)  : %.4e" % pot_hdr)
    print("POT sum of all records                       : %.4e" % rec["pot"].sum())
    fl = rec["flags"]
    for name, m in [("POT cut", 1), ("+ horn", 3), ("+ position", 7), ("+ width (all four)", 15)]:
        sel = (fl & m) == m
        print("POT after %-20s: %.4e   (records %d)" % (name, rec["pot"][sel].sum(), sel.sum()))
    sel = (fl == 15) & inev
    print("POT all four cuts AND event in SPINE list    : %.4e" % rec["pot"][sel].sum())
    if good_runs is not None:
        good = np.array(good_runs)
        grec = np.isin(rec["run"], good); gev = np.isin(ev["run"], good)
        runs = np.unique(ev["run"])
        print("\n--- good runs list: %d runs, %d present in these files, %d files-runs not in list ---" %
              (len(good), len(np.intersect1d(good, runs)), len(np.setdiff1d(runs, good))))
        print("events in good runs: %d of %d" % (gev.sum(), len(gev)))
        print("POT uncut, good runs only                    : %.4e" % rec["pot"][grec].sum())
        sel = (fl == 15) & inev & grec
        print("POT good runs + four cuts + event in list    : %.4e" % rec["pot"][sel].sum())
        print("events passing good runs + all four cuts     : %d" % ((ev["flags"] == 15) & gev).sum())
        print("fraction of uncut POT kept: %.4f" % (rec["pot"][sel].sum() / pot_hdr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir"); ap.add_argument("--filelist"); ap.add_argument("--files", nargs="*")
    ap.add_argument("--max-files", type=int); ap.add_argument("--nproc", type=int, default=4)
    ap.add_argument("--no-width", action="store_true", help="skip the (slow) beam width fit")
    ap.add_argument("--good-runs", help="text file with one good run number per line")
    ap.add_argument("--fit-init", default="howard", choices=["howard", "data"])
    ap.add_argument("--out", default="qc_pilot.npz", help="single output file (non-chunked mode)")
    ap.add_argument("--chunk-dir", help="resumable mode: write one npz per block of files into this directory")
    ap.add_argument("--chunk-size", type=int, default=200)
    ap.add_argument("--merge", help="only merge+summarise the chunk npz files of this directory")
    o = ap.parse_args()
    good = None
    if o.good_runs:
        good = [int(l) for l in open(o.good_runs) if l.strip() and not l.startswith("#")]

    if o.merge:
        cf = sorted(glob.glob(os.path.join(o.merge, "chunk_*.npz")))
        if not cf: sys.exit("no chunk_*.npz files in " + o.merge)
        parts = [load_npz(c) for c in cf]
        ev = {k: np.concatenate([p[0][k] for p in parts]) for k in parts[0][0]}
        rec = {k: np.concatenate([p[1][k] for p in parts]) for k in parts[0][1]}
        secs = np.concatenate([p[3] for p in parts]); bad = sum([p[4] for p in parts], [])
        print("chunks merged:", len(cf))
        summarize(ev, rec, sum(p[2] for p in parts), secs, bad, good, nfiles=len(secs))
        np.savez(os.path.join(o.merge, "merged_tables.npz"), **{"ev_" + k: v for k, v in ev.items()},
                 **{"rec_" + k: v for k, v in rec.items()}, pot_hdr=sum(p[2] for p in parts))
        print("saved", os.path.join(o.merge, "merged_tables.npz"))
        return

    files = o.files or []
    if o.dir: files += sorted(glob.glob(os.path.join(o.dir, "*.flat.root")))
    if o.filelist: files += [l.strip() for l in open(o.filelist) if l.strip()]
    if o.max_files: files = files[:o.max_files]
    if not files: sys.exit("no input files given")
    print("files:", len(files), flush=True)

    if o.chunk_dir:
        os.makedirs(o.chunk_dir, exist_ok=True)
        nblk = (len(files) + o.chunk_size - 1) // o.chunk_size
        t0 = time.time()
        for i in range(nblk):
            out = os.path.join(o.chunk_dir, "chunk_%04d.npz" % i)
            if os.path.exists(out):
                print("block %d/%d exists, skipping" % (i + 1, nblk), flush=True); continue
            r = run_block(files[i * o.chunk_size:(i + 1) * o.chunk_size], o)
            if r is None:
                print("block %d/%d: no readable files, not saved" % (i + 1, nblk), flush=True); continue
            ev, rec, pot, secs, bad = r
            save_npz(out + ".tmp.npz", ev, rec, pot, secs, bad); os.replace(out + ".tmp.npz", out)
            print("block %d/%d done | %.0f min elapsed | block POT %.3e | failed files %d" %
                  (i + 1, nblk, (time.time() - t0) / 60., pot, len(bad)), flush=True)
        os.execvp(sys.executable, [sys.executable, os.path.abspath(__file__), "--merge", o.chunk_dir] +
                  (["--good-runs", o.good_runs] if o.good_runs else []))

    r = run_block(files, o)
    if r is None: sys.exit("no readable files")
    ev, rec, pot, secs, bad = r
    save_npz(o.out, ev, rec, pot, secs, bad)
    summarize(ev, rec, pot, secs, bad, good, nfiles=len(secs))
    print("saved", o.out)


if __name__ == "__main__":
    main()
