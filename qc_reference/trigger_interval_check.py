#!/usr/bin/env python3
"""Gate-versus-spill interval check (Dihan's trigger quality check).

NOTE: the original check was run once as an inline command and was never saved
as a file. The per-run body (d_end / d_start, boundary shift of 100 ms) below is
copied from that command. The database read at the top of the loop and the
final merge are reconstructed. After running, the printed totals must reproduce
the original run: 1,029,152 intervals, 1,024,626 matching (99.560%), 4,526
mismatching, and 930,046 of 930,047 events with a gate row.

For consecutive type-2 majority gates k, k+1 of a run (sorted by event_no):
  diff[k] = (gate_id_numi[k+1] - gate_id_numi[k]) - (number of IFBeam spills in
            (t[k] - 100 ms, t[k+1] - 100 ms])      ; 0 means the interval matches
  d_end[k+1] = diff[k]   interval ENDING at the event   (previous gate -> this gate)
  d_start[k] = diff[k]   interval STARTING at the event (this gate -> next gate)
Gate time (ms, UTC) = (wr_seconds - 37) * 1000 + wr_nanoseconds / 1e6.
Inputs (in the working directory): full_chunks/merged_tables.npz,
ifbeam_spills/spills_run<run>.csv (first column = clock in ms UTC).
Outputs: trigger_interval.npz (per gate) and ev_trigger_interval.npz (per event
of merged_tables.npz, 99 where the event has no gate row).
"""
import collections, getpass
import numpy as np
import pg8000.dbapi as pg

T = np.load("full_chunks/merged_tables.npz")
runs = sorted(int(r) for r in np.unique(T["ev_run"]))
conn = pg.connect(host="ifdbdaqrep01.fnal.gov", port=5455, database="icarus_trigger_prd",
                  user="triggerdb_reader", password=getpass.getpass("trigger DB password: "))
cur = conn.cursor()
Rl, El, Dend, Dstart = [], [], [], []
hist = collections.Counter(); perrun = []; nonmono = 0
shifts = {30.0: 0, 100.0: 0, 500.0: 0, 1000.0: 0}
for run in runs:
    cur.execute("SELECT event_no, wr_seconds, wr_nanoseconds, gate_id_numi FROM triggerdata "
                "WHERE run_number=%s AND gate_type=2 AND trigger_type=0 ORDER BY event_no", (run,))
    rows = cur.fetchall()
    ev = np.array([r[0] for r in rows], dtype=np.int64)
    t = np.array([(r[1] - 37) * 1000.0 + r[2] / 1e6 for r in rows])
    numi = np.array([r[3] for r in rows], dtype=np.int64)
    s = np.sort(np.loadtxt("ifbeam_spills/spills_run%d.csv" % run, delimiter=",", skiprows=1, ndmin=2)[:, 0])
    # ---- from here on: copied from the original command ----
    if not np.all(np.diff(t) > 0): nonmono += 1
    dn = numi[1:] - numi[:-1]
    for sh in shifts:
        b = t - sh
        shifts[sh] += int(((dn - (np.searchsorted(s, b[1:], side="right") - np.searchsorted(s, b[:-1], side="right"))) != 0).sum())
    b = t - 100.0
    diff = dn - (np.searchsorted(s, b[1:], side="right") - np.searchsorted(s, b[:-1], side="right"))
    de = np.zeros(len(t), dtype=np.int16); ds = np.zeros(len(t), dtype=np.int16); de[1:] = diff; ds[:-1] = diff
    Rl.append(np.full(len(t), run, dtype=np.int64)); El.append(ev); Dend.append(de); Dstart.append(ds)
    hist.update(diff.tolist()); perrun.append((int((diff != 0).sum()), run, len(diff)))
conn.close()
R = np.concatenate(Rl); E = np.concatenate(El); De = np.concatenate(Dend); Ds = np.concatenate(Dstart)
np.savez("trigger_interval.npz", run=R, event=E, d_end=De, d_start=Ds)
ni = sum(hist.values())
print("intervals:", ni, "| gates == spills: %d (%.3f%%)" % (hist[0], 100.0 * hist[0] / ni))
print("(gates - spills) -> count:", sorted(hist.items(), key=lambda kv: -kv[1])[:12])
print("mismatching intervals for boundary shift 30/100/500/1000 ms:", [shifts[k] for k in sorted(shifts)], "| runs with non-monotonic times:", nonmono)
# ---- merge onto the events of merged_tables.npz (99 = no gate row) ----
kg = (R << 40) | E; ke = (T["ev_run"].astype(np.int64) << 40) | T["ev_evt"].astype(np.int64)
o = np.argsort(kg); i = np.clip(np.searchsorted(kg[o], ke), 0, len(kg) - 1); has = kg[o][i] == ke
d_end = np.where(has, De[o][i], 99).astype(np.int16); d_start = np.where(has, Ds[o][i], 99).astype(np.int16)
np.savez("ev_trigger_interval.npz", run=T["ev_run"], evt=T["ev_evt"], d_end=d_end, d_start=d_start)
print("our events with a gate row:", int(has.sum()), "of", len(ke))
