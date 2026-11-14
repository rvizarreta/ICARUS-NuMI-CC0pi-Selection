import argparse, os, re, ssl, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import numpy as np

TABLES = "full_chunks/merged_tables.npz"
GATES = "trigger_gates.npz"
SPDIR = "ifbeam_spills"
FAILED = "ifbeam_failed_windows.txt"
BASE = "https://dbdata3vm.fnal.gov:9443/ifbeam/data//data?"
STEP = 600000.0
CTX = ssl._create_unverified_context()
RX = re.compile(r'clock="(\d+)">\s*([-+0-9.eE]+)')

def runs_in_tables():
    return sorted(int(r) for r in np.unique(np.load(TABLES)["ev_run"]))

def do_gates():
    import getpass
    import pg8000.dbapi as pg
    pw = getpass.getpass("trigger DB password: ")
    conn = pg.connect(host="ifdbdaqrep01.fnal.gov", port=5455, database="icarus_trigger_prd",
                      user="triggerdb_reader", password=pw)
    cur = conn.cursor()
    R, E, T = [], [], []
    runs = runs_in_tables()
    for n, run in enumerate(runs, 1):
        cur.execute("SELECT event_no, wr_seconds, wr_nanoseconds FROM triggerdata "
                    "WHERE run_number=%s AND gate_type=2 AND trigger_type=0 ORDER BY event_no", (run,))
        for e, s, ns in cur.fetchall():
            R.append(run); E.append(e); T.append((s - 37) * 1000.0 + ns / 1e6)
        if n % 10 == 0: print("runs read:", n, "of", len(runs), "| gates:", len(R), flush=True)
    conn.close()
    np.savez(GATES, run=np.array(R, dtype=np.int64), event=np.array(E, dtype=np.int64), t_ms=np.array(T))
    print("saved", GATES, "with", len(R), "gates")

def stamp(ms):
    from zoneinfo import ZoneInfo
    s = datetime.fromtimestamp(ms / 1000.0, ZoneInfo("America/Chicago")).strftime("%Y-%m-%dT%H:%M:%S.%f%z")
    return s[:-2] + ":" + s[-2:]

def get(t0, t1):
    url = BASE + "v=E:TRTGTD&e=e,a9&t0=%s&t1=%s&f=xml" % (stamp(t0), stamp(t1))
    for k in range(5):
        try:
            body = urllib.request.urlopen(url, timeout=60, context=CTX).read().decode()
            out = [(int(c), float(v)) for c, v in RX.findall(body)]
            for c, _ in out:
                if c < t0 - 2000 or c > t1 + 2000:
                    sys.exit("IFBeam returned a spill outside the requested window: time zone problem")
            time.sleep(0.2)
            return out
        except SystemExit:
            raise
        except Exception:
            time.sleep(3 * (k + 1))
    return None

def do_fetch():
    try: stamp(0)
    except Exception: sys.exit("no time zone data available (zoneinfo)")
    g = np.load(GATES); grun = g["run"]; gt = g["t_ms"]
    os.makedirs(SPDIR, exist_ok=True); t00 = time.time()
    runs = [int(r) for r in np.unique(grun)]
    with ThreadPoolExecutor(4) as ex:
        for n, run in enumerate(runs, 1):
            fn = "%s/spills_run%d.csv" % (SPDIR, run)
            if os.path.exists(fn): continue
            t = gt[grun == run]; lo, hi = float(t.min()) - 5000.0, float(t.max()) + 5000.0
            wins = []; a = lo
            while a < hi: wins.append((a, min(a + STEP, hi))); a += STEP
            res = list(ex.map(lambda w: get(w[0], w[1]), wins))
            spills = {}; nfail = 0
            for w, r in zip(wins, res):
                if r is None:
                    nfail += 1
                    open(FAILED, "a").write("%d %.0f %.0f\n" % (run, w[0], w[1]))
                else:
                    for c, v in r: spills[c] = v
            with open(fn + ".tmp", "w") as f:
                f.write("clock_ms,trtgtd_E12\n")
                for c in sorted(spills): f.write("%d,%.4f\n" % (c, spills[c]))
            os.replace(fn + ".tmp", fn)
            print("run %d: %d windows, %d spills, %d failed windows | run %d of %d | %.0f s"
                  % (run, len(wins), len(spills), nfail, n, len(runs), time.time() - t00), flush=True)

def do_analyze(thr):
    g = np.load(GATES); grun = g["run"]; gev = g["event"]; gt = g["t_ms"]
    d = np.full(len(gt), np.inf)
    for run in [int(r) for r in np.unique(grun)]:
        fn = "%s/spills_run%d.csv" % (SPDIR, run)
        if not os.path.exists(fn): print("no spill file for run", run); continue
        s = np.sort(np.loadtxt(fn, delimiter=",", skiprows=1, ndmin=2)[:, 0])
        if len(s) < 2: continue
        m = np.where(grun == run)[0]; t = gt[m]
        i = np.clip(np.searchsorted(s, t), 1, len(s) - 1)
        near = np.where(np.abs(t - s[i - 1]) <= np.abs(t - s[i]), s[i - 1], s[i])
        d[m] = t - near
    np.savez("trigger_dt.npz", run=grun, event=gev, d_ms=d)
    bad = np.abs(d) > thr
    print("gates:", len(d), "| |gate - spill| > %g ms: %d (%.4f%%)" % (thr, bad.sum(), 100.0 * bad.mean()))
    print("percentiles 1,50,99 of gate - spill for |d| <= thr (ms):", np.round(np.percentile(d[~bad], [1, 50, 99]), 2))
    runs = np.unique(grun)
    cnt = sorted(((int(bad[grun == r].sum()), int(r), int((grun == r).sum())) for r in runs), reverse=True)
    print("10 runs with most failing gates (n_fail, run, n_gates):", cnt[:10])
    T = np.load(TABLES); k_ours = (T["ev_run"].astype(np.int64) << 40) | T["ev_evt"].astype(np.int64)
    k_bad = (grun[bad].astype(np.int64) << 40) | gev[bad].astype(np.int64)
    fo = np.isin(k_ours, k_bad); ok15 = T["ev_flags"] == 15
    print("our events:", len(k_ours), "| failing the timing check:", int(fo.sum()), "(%.4f%%)" % (100.0 * fo.mean()),
          "| of those, passing the beam cuts:", int((fo & ok15).sum()))
    if os.path.exists(FAILED): print("IFBeam windows that failed after retries:", sum(1 for _ in open(FAILED)))

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gates", action="store_true"); ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--analyze", action="store_true"); ap.add_argument("--thr-ms", type=float, default=100.0)
    a = ap.parse_args()
    if a.gates: do_gates()
    if a.fetch: do_fetch()
    if a.analyze: do_analyze(a.thr_ms)
