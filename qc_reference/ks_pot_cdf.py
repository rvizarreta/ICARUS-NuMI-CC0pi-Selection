#!/usr/bin/env python3
"""
KS-style test: accumulated POT vs time against accumulated selected candidates vs time
(Howard et al. Sec. 8.2.2 procedure).

Inputs
  --tables   merged_tables.npz from make_qc_lists.py (needs rec_time and ev_time, see patch)
  --cands    npz with arrays run, subrun, evt for the selected candidates (one file per selection)
  --good-runs  optional text file, one run per line

POT CDF   : spills passing all four QC cuts (flags==15), in good runs, whose triggering event
            is in the SPINE event list; ordered in time, cumulative POT normalised to 1.
Cand. CDF : candidate times sorted, i/N.
D         : max over candidates of |candCDF(t_i) - potCDF(t_i)|, potCDF evaluated at the spill
            closest in time to the candidate.
p-value   : kstwo survival function at D with N_eff = Ns*Nc/(Ns+Nc) (rounded to an integer,
            scipy's kstwo needs integer n), plus a pseudo-experiment cross-check.
"""
import argparse
import numpy as np
from scipy.stats import kstwo


def key(run, sub, evt):
    return (run.astype(np.int64) << 40) | (sub.astype(np.int64) << 26) | evt.astype(np.int64)


def load_tables(path, good_runs=None, qc="all"):
    d = np.load(path)
    ev = {k[3:]: d[k] for k in d.files if k.startswith("ev_")}
    rec = {k[4:]: d[k] for k in d.files if k.startswith("rec_")}
    for need, tab in (("time", ev), ("time", rec)):
        if need not in tab:
            raise SystemExit("tables have no 'time' column: add spill time to make_qc_lists.py first")
    kev = key(ev["run"], ev["subrun"], ev["evt"])
    krec = key(rec["run"], rec["subrun"], rec["event"])
    if qc == "all":      # all four QC cuts and event in the SPINE list
        sel = (rec["flags"] == 15) & np.isin(krec, kev)
    else:                # qc == "none": every spill record, no quality cuts (first pass)
        sel = np.ones(len(krec), bool)
    if good_runs is not None:
        sel &= np.isin(rec["run"], good_runs)
    t, pot, run = rec["time"][sel], rec["pot"][sel], rec["run"][sel]
    o = np.argsort(t)
    return ev, t[o], pot[o], run[o]


def candidate_times(ev, cands):
    """look the candidate (run, subrun, evt) up in the event table to get the trigger-spill time"""
    kev = key(ev["run"], ev["subrun"], ev["evt"])
    kc = key(cands["run"], cands["subrun"], cands["evt"])
    idx = {k: i for i, k in enumerate(kev)}
    miss = [k for k in kc if k not in idx]
    if miss:
        print("WARNING: %d candidates not found in event table (dropped)" % len(miss))
    return np.sort(np.array([ev["time"][idx[k]] for k in kc if k in idx]))


def ks_test(t_spill, pot, t_cand, n_toys=1000, seed=1):
    pot_cdf = np.cumsum(pot) / pot.sum()
    n_c, n_s = len(t_cand), len(t_spill)
    # nearest spill to each candidate
    j = np.clip(np.searchsorted(t_spill, t_cand), 1, n_s - 1)
    j -= (t_cand - t_spill[j - 1]) < (t_spill[j] - t_cand)
    cand_cdf = np.arange(1, n_c + 1) / n_c
    diff = np.abs(cand_cdf - pot_cdf[j])
    D = diff.max()
    n_eff = n_s * n_c / (n_s + n_c)
    p = kstwo.sf(D, int(round(n_eff)))
    # pseudo-experiments: draw n_c spills with probability ~ POT
    rng = np.random.default_rng(seed)
    w = pot / pot.sum()
    Dt = np.empty(n_toys)
    for i in range(n_toys):
        k = np.sort(rng.choice(n_s, size=n_c, p=w))
        Dt[i] = np.abs(cand_cdf - pot_cdf[k]).max()
    p_toy = (Dt >= D).mean()
    return dict(D=D, n_spills=n_s, n_cand=n_c, n_eff=n_eff, p_scipy=p, p_toys=p_toy,
                dt=t_cand - t_spill[j], diff=diff, pot_cdf=pot_cdf, cand_cdf=cand_cdf,
                t_spill=t_spill, t_cand=t_cand)


def plots(res, label, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    ax[0].plot(res["t_spill"], res["pot_cdf"], label="POT")
    ax[0].plot(res["t_cand"], res["cand_cdf"], label="candidates")
    ax[0].set_xlabel("time [Unix s]"); ax[0].set_ylabel("CDF"); ax[0].legend(); ax[0].set_title(label)
    ax[1].plot(res["t_cand"], res["diff"]); ax[1].set_xlabel("time [Unix s]")
    ax[1].set_ylabel("|candidate CDF - POT CDF|")
    ax[1].set_title("D = %.4f, p = %.3f" % (res["D"], res["p_scipy"]))
    ax[2].hist(np.clip(res["dt"], -20, 20), bins=80); ax[2].set_yscale("log")
    ax[2].set_xlabel("candidate time - matched spill time [s]"); ax[2].set_title("matching check")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)


def per_run(run_spill, pot, cand_run, label, out):
    """Supplementary check: candidates per POT in each run, compared with the constant average rate.
    Chi2 = sum (N - rate*POT)^2 / (rate*POT) over runs, ndf = nruns - 1 (Poisson, approximate when
    the expected count in a run is small)."""
    from scipy.stats import chi2 as chi2d
    runs = np.unique(run_spill)
    pot_r = np.array([pot[run_spill == r].sum() for r in runs])
    n_r = np.array([(cand_run == r).sum() for r in runs])
    n_out = int((~np.isin(cand_run, runs)).sum())
    keep = pot_r > 0
    runs, pot_r, n_r = runs[keep], pot_r[keep], n_r[keep]
    rate = n_r.sum() / pot_r.sum()
    exp = rate * pot_r
    chi2 = ((n_r - exp) ** 2 / exp).sum(); ndf = len(runs) - 1
    p = chi2d.sf(chi2, ndf)
    print("   per-run: %d runs, mean rate %.3e cand/POT, chi2/ndf = %.1f/%d, p = %.3f "
          "(runs with expected <5 events: %d; candidates in runs with no POT: %d)"
          % (len(runs), rate, chi2, ndf, p, (exp < 5).sum(), n_out))
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.errorbar(runs, n_r / pot_r / rate, yerr=np.sqrt(np.maximum(n_r, 1)) / pot_r / rate, fmt="o", ms=3)
    ax.axhline(1, color="k", lw=1)
    ax.set_xlabel("run"); ax.set_ylabel("(candidates / POT) / mean")
    ax.set_title("%s: chi2/ndf = %.1f/%d, p = %.3f" % (label, chi2, ndf, p))
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return dict(chi2=chi2, ndf=ndf, p=p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", required=True)
    ap.add_argument("--cands", nargs="+", required=True, help="name=file.npz, one per selection")
    ap.add_argument("--good-runs"); ap.add_argument("--toys", type=int, default=1000)
    ap.add_argument("--qc", choices=["all", "none"], default="all",
                    help="all = spills passing the 4 QC cuts (final); none = uncut POT (first pass)")
    ap.add_argument("--tag", default="", help="suffix for output figure names")
    o = ap.parse_args()
    good = None
    if o.good_runs:
        good = np.array([int(l) for l in open(o.good_runs) if l.strip() and not l.startswith("#")])
    ev, t_spill, pot, run_spill = load_tables(o.tables, good, o.qc)
    print("POT CDF (qc=%s): %d spills, %.4e POT" % (o.qc, len(t_spill), pot.sum()))
    for spec in o.cands:
        name, path = spec.split("=")
        c = np.load(path)
        res = ks_test(t_spill, pot, candidate_times(ev, c), o.toys)
        print("%-12s Ncand=%d  Neff=%.1f  D=%.4f  p(kstwo)=%.3f  p(toys)=%.3f" %
              (name, res["n_cand"], res["n_eff"], res["D"], res["p_scipy"], res["p_toys"]))
        plots(res, name, "ks_%s%s.png" % (name, o.tag))
        per_run(run_spill, pot, c["run"], name, "perrun_%s%s.png" % (name, o.tag))


if __name__ == "__main__":
    main()
