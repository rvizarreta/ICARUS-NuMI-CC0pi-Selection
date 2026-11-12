#!/usr/bin/env python3
"""Look at real MTGTDS beam-profile channels and compare the two fit start-value choices.
usage: python3 width_check.py FILE.flat.root [N]"""
import sys, numpy as np, uproot, awkward as ak
import make_qc_lists as q
f = uproot.open(sys.argv[1]); t = f[[k for k in f.keys() if 'recTree' in k or 'Events' in k][0]]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
a = t.arrays(['rec.hdr.numiinfo.MTGTDS', 'rec.hdr.numiinfo.MTGTDS..length'], entry_stop=1)
flat = ak.to_numpy(ak.flatten(a['rec.hdr.numiinfo.MTGTDS'], axis=None))
lens = ak.to_numpy(ak.flatten(a['rec.hdr.numiinfo.MTGTDS..length'], axis=None))
m, ok = q.to2d(flat, lens, q.NMTGT)
print('records', len(m), '| channel counts', set(lens.tolist()))
c = -m[0, 103:151]
print('record 0, h_chan = -MTGTDS[103:151]: min %.3g max %.3g mean %.3g' % (c.min(), c.max(), c.mean()))
print('h_chan:', np.round(c, 2).tolist())
res = {'howard': [], 'data': []}
for mode in res:
    for i in range(min(N, len(m))):
        res[mode].append((q.sigma_fit(-m[i, 103:151], mode), q.sigma_fit(-m[i, 151:199], mode)))
    r = np.array(res[mode])
    inr = ((r > 0.57) & (r <= 1.88)).all(1)
    print('%-6s start values: sigma_H median %.3f, sigma_V median %.3f | in (0.57,1.88]: %d of %d | NaN: %d' %
          (mode, np.nanmedian(r[:, 0]), np.nanmedian(r[:, 1]), inr.sum(), len(r), np.isnan(r).any(1).sum()))
