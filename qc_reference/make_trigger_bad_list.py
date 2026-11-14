#!/usr/bin/env python3
"""Write the bad-event list read by trigger_quality_cut (event_cuts.h).

An event is bad if its interval check fails. d_end / d_start come from
ev_trigger_interval.npz (made by trigger_qc.py on the GPVM). For consecutive
majority gates k and k+1 of a run, diff[k] = (change of gate_id_numi) - (number
of IFBeam spills in the interval); 0 means the interval matches. trigger_qc.py
stores it as d_end[k+1] = diff[k] (the interval ENDING at the event: previous
gate -> this gate) and d_start[k] = diff[k] (the interval STARTING at the
event: this gate -> next gate). An event with no gate row has 99 in both.

usage: make_trigger_bad_list.py <ev_trigger_interval.npz> <out.txt> --policy {end,start,either}
"""
import argparse
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("npz")
ap.add_argument("out")
ap.add_argument("--policy", required=True, choices=["end", "start", "either"])
a = ap.parse_args()

t = np.load(a.npz)
bad_end = t["d_end"] != 0
bad_start = t["d_start"] != 0
bad = {"end": bad_end, "start": bad_start, "either": bad_end | bad_start}[a.policy]
run = t["run"].astype(np.int64)
evt = t["evt"].astype(np.int64)
with open(a.out, "w") as f:
    for r, e in zip(run[bad], evt[bad]):
        f.write(f"{r} {e}\n")
print(f"events: {len(run)}  bad ({a.policy}): {bad.sum()}  written to {a.out}")
