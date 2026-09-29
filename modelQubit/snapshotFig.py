# -*- coding: utf-8 -*-
"""
snapshotFig.py -- three-panel snapshot figure of one trial for the paper
(replaces results/phase1-3.png). Uses mqAnalysis.py's loader and trial table.

  python snapshotFig.py ../modelQubit_09-24-2026_220244.log 46 trialSnapshots.pdf

Panels: (a) force off, state resting at the previous measurement's pole;
(b) new axis set, force on, state starting to move; (c) final equilibrium.
Coordinates are logged x, y divided by the median radius (natural units).
"""
import math
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

import mqAnalysis as M

LOG = sys.argv[1] if len(sys.argv) > 1 else '../modelQubit_09-24-2026_220244.log'
K = int(sys.argv[2]) if len(sys.argv) > 2 else 46
OUT = sys.argv[3] if len(sys.argv) > 3 else 'trialSnapshots.pdf'
N_MID = 2          # samples after force-on shown in panel (b)
TRAIL = 12         # samples of trail drawn behind the state

d = M.load_log(LOG)
trials, _ = M.segment(d)
T = M.trial_table(d, trials)
tr = trials[K]
R = float(np.median(d['r']))


def axis_at(i):
    return math.atan2(d['s'][i], d['c'][i])


# previous measurement axis (the one the state collapsed onto last time)
prevAxis = axis_at(trials[K - 1]['onStart']) if K > 0 else None
newAxis = axis_at(tr['onStart'])
iPre, iOn, iEnd = tr['pre'], tr['onStart'], tr['onEnd']
frames = [(iPre, prevAxis, False, '(a) Initial state'),
          (min(iOn + N_MID, iEnd), newAxis, True, '(b) New axis, force on'),
          (iEnd, newAxis, True, '(c) Final state')]

fig, axs = plt.subplots(1, 3, figsize=(3.5, 1.45))
for ax, (i, axisAng, dashedToo, title) in zip(axs, frames):
    ax.set_aspect('equal')
    ax.add_patch(Circle((0, 0), 1, fill=False, ec='0.6'))
    if axisAng is not None:
        ca, sa = math.cos(axisAng), math.sin(axisAng)
        ax.plot([0, ca], [0, sa], color='tab:orange', lw=1.3)
        ax.plot([0, -ca], [0, -sa], color='tab:orange', lw=1.3, ls='--')
        ax.plot(ca, sa, 'o', color='tab:orange', mfc='white', ms=5)
        ax.text(1.3 * ca, 1.3 * sa, r'$\uparrow$', ha='center', va='center', fontsize=9)
        ax.text(-1.3 * ca, -1.3 * sa, r'$\downarrow$', ha='center', va='center', fontsize=9)
    lo = max(0, i - TRAIL)
    ax.plot(d['x'][lo:i + 1] / R, d['y'][lo:i + 1] / R, '-', color='tab:blue',
            alpha=0.35, lw=1.3)
    ax.plot(d['x'][i] / R, d['y'][i] / R, 'o', color='tab:blue', ms=5)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=7)
fig.tight_layout(pad=0.3, w_pad=0.2)
fig.savefig(OUT)
print(f'trial {K}: D = {T["D_deg"][K]:.1f} deg, P_Born = {T["P_born"][K]:.3f}, '
      f'w = {T["w"][K]:.3f}, outcome = {T["outcome"][K]}, t(on) = {d["t"][iOn]:.2f} s, '
      f'panel (b) at +{d["t"][frames[1][0]] - d["t"][iOn]:.2f} s, '
      f'panel (c) at +{d["t"][iEnd] - d["t"][iOn]:.2f} s')
