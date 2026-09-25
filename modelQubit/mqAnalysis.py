# -*- coding: utf-8 -*-
"""
mqAnalysis.py -- Model Qubit data analysis
Benjamin Krawciw

Replaces mqDataAnalysis.py. One log in, three questions out:

  GOAL 1  Does the circuit faithfully implement the ODE model?
          1a  theta transfer:   logged theta vs the designed theta(u, phi)
          1b  radius:           settled |r| vs the designed sqrt(14/19)
          1c  dynamics:         first-order tangential law r dpsi/dt = g theta (n.t)
                                (g vs the schematic's 50 /s; drift F0; forward prediction)
          1d  sign rule:        final pole == sign(theta_0) * (c, s)
          1e  end to end:       outcome vs the ideal model's outcome for the same
                                (D, phi); threshold-noise fit (mu, sigma_m)

  GOAL 2  Do the empirical probabilities obey the Born rule?
          Per-bin error budget built from the data (see ERROR MODEL below),
          chi^2 and p-value against cos^2(D/2) evaluated at the realized angles.

  GOAL 3  Animations: one per trial (diagnostic panels) and an overview
          animation of every trial collapsing at once (slides).

Outputs go to OUTDIR (default: analysis_<log name>/ next to this script):
  summary.txt, trials.csv, trialCycle.pdf, empiricalProb.pdf, odeChecks.pdf,
  animations/trial_###.gif (or .mp4), animations/overview.gif

Runs from Spyder (edit CONFIG) or the command line:
  python mqAnalysis.py ../modelQubit_09-24-2026_200144.log --anim failures

-------------------------------------------------------------------------------
ERROR MODEL FOR THE EMPIRICAL PROBABILITIES
-------------------------------------------------------------------------------
Trial i has realized relative angle D_i (axis vs the state just before force-on)
and hidden variable phi_i, normalized to w_i in [0, 1] by the designed sweep.
The ideal model says  up  <=>  w_i < P_i := cos^2(D_i / 2).

For bin j with n_j trials, p_hat_j = (ups + 0.5 * ambiguous) / n_j is compared
with the Born value at the realized angles, P_j = mean_i P_i (so any difference
between commanded and realized angle is accounted for, not added to the bar).

Error bar = spread of p_hat_j under the hypothesis "this circuit, this protocol".
Two ingredients:
  * Device threshold noise. Fit P(up_i) = Phi((m_i - mu)/sigma_m) by maximum
    likelihood over all trials, where m_i = P_i - w_i is the margin from the
    ideal threshold. sigma_m is the circuit's threshold jitter and mu its bias,
    both in probability units.
  * Hidden-variable sampling. phi takes N_phi discrete values, so p_hat is a
    quadrature of a step function. Treat each phi value as one uniform draw
    w_c + U from its grid cell (width h). Repeat suites reuse the same phi, so
    they share U and add no new hidden-variable samples.
Given U, trial i in cell c is up with probability p_i(U) = Phi((P_i - w_c - U - mu)
/ sigma_m). The law of total variance splits each cell's up-count variance into
  sigma_HV^2  = sum_c Var_U[ sum_i p_i(U) ]          / n_j^2   (hidden variable)
  sigma_dev^2 = sum_c E_U [ sum_i p_i(1 - p_i) ]      / n_j^2   (circuit noise)
plus sigma_amb^2 = 0.25 n_amb / n_j^2 for ambiguous outcomes counted as 1/2.
sigma_j^2 = sigma_HV^2 + sigma_dev^2 + sigma_amb^2, floored at SIGMA_FLOOR (only binds
where the phi grid is symmetric about the threshold, e.g. D = 90 deg).
A Monte Carlo of this model reproduces the spread of p_hat to within a few
percent in every bin (checked on the 09-24-2026 log).

Also reported:
  ideal     E_j = mean_i [w_i < P_i]: a perfect circuit on exactly these phi values.
  gridBias  E_U[p_hat] - P_j for a perfect circuit: the quadrature bias of the
            phi grid itself. An endpoint-inclusive grid (phiMin + i*step/(N-1))
            is biased toward 1/2 by up to ~1/(2N); a midpoint grid
            (phiMin + (i + 0.5)*step/N) is unbiased.
"""

import argparse
import csv
import math
import os
import re

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation as ani
from matplotlib.patches import Circle


# =============================================================================
# CONFIG
# =============================================================================
LOGFILE = '../modelQubit_09-24-2026_220244.log'
OUTDIR = None           # None -> analysis_<log stem>/ beside this script
ANIMATE = 'failures'    # 'none' | 'failures' | 'all' | list of trial numbers
ANIM_FORMAT = 'gif'     # 'gif' or 'mp4' (mp4 needs ffmpeg; better for slides)
ANIM_SPEED = 1.0        # 1.0 = real time
OVERVIEW = True         # overview animation of all trials at once
OVERWRITE = False       # re-render animations that already exist
TRIALCYCLE = (3, 5)     # first/last trial shown in trialCycle.pdf
N_BOOT = 400            # bootstrap resamples (by trial) for parameter errors
SEED = 0

# Circuit design (must match modelQubit.ino and the paper's circuit analysis)
THETA_GAIN = 10.0 / 47.0                  # (50/47 V) / 5 V
R_EQ = math.sqrt(14.0 / 19.0)             # radius set by the 1 - r^2 circuit
PHI_MIN = THETA_GAIN * (1.0 - R_EQ)       # designed sweep, fractions of 5 V
PHI_MAX = THETA_GAIN * (1.0 + R_EQ)
G_DESIGN = 5.0 / (10.0 * 10e-3)           # AD633 /10 V, tau = 10 ms -> 50 /s
N_ANGLES = 10                             # firmware m
ANGLE_MAX = np.pi / 2
C_OFF = 2.0 * 128.0 / 255.0 - 1.0         # measureOff() writes PWM 128

# Tolerances from component tolerances (toleranceMC.py: resistors +/-1%, capacitors
# +/-10%, AD633 scale +/-0.25%, propagated through the circuit diagrams). Each entry is
# the central 95% interval of (value/design - 1), or of mu in probability units:
#   'typ' = AD633 offsets at datasheet typical (5 mV)   -> outside: WARN
#   'max' = AD633 offsets at datasheet maximum (50/30 mV) -> outside: FAIL
# Not included: the +5 V rail. b and g scale as V_rail, R as V_rail^-1/2 (see toleranceMC.py).
TOL = {
    'a0':   {'typ': (-0.036, 0.036), 'max': (-0.048, 0.049)},
    'b':    {'typ': (-0.018, 0.018), 'max': (-0.018, 0.018)},
    'gphi': {'typ': (-0.021, 0.021), 'max': (-0.021, 0.022)},
    'R':    {'typ': (-0.010, 0.010), 'max': (-0.024, 0.024)},
    'g':    {'typ': (-0.092, 0.110), 'max': (-0.092, 0.111)},   # position-integrator C
    'mu':   {'typ': (-0.011, 0.013), 'max': (-0.024, 0.030)},   # Born threshold bias
}
P_VALUE_MIN = 0.05         # Born chi^2 p-value must exceed this
SIGMA_FLOOR = 0.002        # minimum per-bin error bar (well below the ~0.01 bars of a dense grid)
OUTCOME_TOL_DEG = 20.0     # final state within this of a pole is a definite outcome;
                           # otherwise the collapse was incomplete -> ambiguous (1/2)


# =============================================================================
# SMALL MATH HELPERS (numpy only; no scipy needed)
# =============================================================================
def wrap(a):
    """Wrap angle(s) to (-pi, pi]."""
    return (np.asarray(a) + np.pi) % (2.0 * np.pi) - np.pi


def normcdf(z):
    """Standard normal CDF (Abramowitz-Stegun 7.1.26, |err| < 1.5e-7)."""
    z = np.asarray(z, dtype=float)
    x = np.abs(z) / np.sqrt(2.0)
    tt = 1.0 / (1.0 + 0.3275911 * x)
    poly = tt * (0.254829592 + tt * (-0.284496736 + tt * (1.421413741
            + tt * (-1.453152027 + tt * 1.061405429))))
    erf = 1.0 - poly * np.exp(-x * x)
    return 0.5 * (1.0 + np.sign(z) * erf)


def chi2_sf(x, k):
    """P(chi^2_k > x): regularized upper incomplete gamma Q(k/2, x/2)."""
    a, xx = 0.5 * k, 0.5 * x
    if xx <= 0:
        return 1.0
    if xx < a + 1.0:                       # series for P, return 1 - P
        term = total = 1.0 / a
        n = a
        for _ in range(500):
            n += 1.0
            term *= xx / n
            total += term
            if abs(term) < abs(total) * 1e-14:
                break
        return 1.0 - total * math.exp(-xx + a * math.log(xx) - math.lgamma(a))
    b = xx + 1.0 - a                       # continued fraction for Q
    c = 1.0 / 1e-300
    d = 1.0 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1.0 / d
        h *= d * c
        if abs(d * c - 1.0) < 1e-14:
            break
    return math.exp(-xx + a * math.log(xx) - math.lgamma(a)) * h


def rsq(y, yhat):
    ss = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - np.sum((y - yhat) ** 2) / ss if ss > 0 else float('nan')


# =============================================================================
# LOAD AND SEGMENT
# =============================================================================
def load_log(fname):
    """Read a serial log. Tolerates terminal-logger banners (e.g. PuTTY) and partial
    lines before the firmware header: the header is the line containing 'Voltages
    are recorded', the column line follows it, and only rows of exactly 8 numbers
    after that are kept."""
    with open(fname, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.read().splitlines()
    hIdx = next((i for i, ln in enumerate(lines) if 'Voltages are recorded' in ln), None)
    if hIdx is None:
        raise ValueError(f'{fname}: no firmware header line found')
    header = lines[hIdx]
    vals = []
    for ln in lines[hIdx + 2:]:
        fields = [v.strip() for v in ln.split(',') if v.strip()]
        if len(fields) != 8:
            continue
        try:
            vals.append([float(v) for v in fields])
        except ValueError:
            continue
    arr = np.array(vals)
    m = re.search(r'fw v(\d+)', header)
    fw = int(m.group(1)) if m else 1
    d = dict(zip(('t', 'c', 's', 'on', 'th', 'ph', 'x', 'y'), arr.T))
    d['on'] = d['on'].astype(bool)
    if fw < 2:
        # pre-v2 firmware logged theta with the x/y conversion; see mqDataAnalysis.py
        d['th'] = 0.625 * d['th'] + 0.125
    d['fw'] = fw
    d['thetaLag'] = 0 if fw >= 2 else 2
    d['psi'] = np.arctan2(d['y'], d['x'])
    d['r'] = np.hypot(d['x'], d['y'])
    d['cEff'] = np.where(d['on'], d['c'], C_OFF)
    d['sEff'] = np.where(d['on'], d['s'], C_OFF)
    return d


def segment(d):
    """Trials start at each force-on rising edge (after the axis is first set)."""
    on, c, s = d['on'], d['c'], d['s']
    armed = (c != 0.0) | (s != 0.0)
    starts = [i for i in range(1, len(on)) if on[i] and not on[i - 1] and armed[i]]
    ends = starts[1:] + [len(on)]
    trials = []
    for a, b in zip(starts, ends):
        onIdx = a + np.flatnonzero(on[a:b])
        # contiguous force-on block at the start of the trial
        brk = np.flatnonzero(np.diff(onIdx) > 1)
        onIdx = onIdx[:brk[0] + 1] if len(brk) else onIdx
        trials.append(dict(start=a, end=b - 1, onStart=onIdx[0], onEnd=onIdx[-1],
                           pre=a - 1 if not on[a - 1] else a))
    trialOf = np.full(len(on), -1)
    for k, tr in enumerate(trials):
        trialOf[tr['start']:tr['end'] + 1] = k
    return trials, trialOf


def trial_table(d, trials):
    """One row per trial with everything needed for goals 1-3."""
    rows = []
    for k, tr in enumerate(trials):
        i0, i1, ip = tr['onStart'], tr['onEnd'], tr['pre']
        axis = math.atan2(d['s'][i0], d['c'][i0])
        D = float(wrap(axis - d['psi'][ip]))          # firmware rotates +D
        th0 = float(d['th'][min(i0 + d['thetaLag'], i1)])
        fin = abs(math.degrees(wrap(d['psi'][i1] - axis)))
        outcome = (1.0 if fin <= OUTCOME_TOL_DEG else
                   0.0 if fin >= 180.0 - OUTCOME_TOL_DEG else 0.5)
        phi = float(d['ph'][i0])
        w = (phi - PHI_MIN) / (PHI_MAX - PHI_MIN)
        P = math.cos(0.5 * D) ** 2
        rows.append(dict(
            trial=k, t0=float(d['t'][i0]), D=D, D_deg=math.degrees(D),
            bin=int(round(D / (ANGLE_MAX / (N_ANGLES - 1)))),
            phi=phi, w=w, P_born=P, margin=P - w, theta0=th0,
            final_deg=math.degrees(wrap(d['psi'][i1] - axis)),
            outcome=outcome, sign_pred=1.0 if th0 >= 0 else 0.0,
            ideal_up=1.0 if w < P else 0.0,
        ))
    T = {key: np.array([r[key] for r in rows]) for key in rows[0]}
    T['sign_ok'] = T['outcome'] == T['sign_pred']
    T['model_ok'] = T['outcome'] == T['ideal_up']
    T['ambiguous'] = T['outcome'] == 0.5
    return T


# =============================================================================
# GOAL 1: ODE FIDELITY
# =============================================================================
def boot_by_trial(fn, groups, nBoot, rng):
    """Bootstrap a statistic by resampling whole trials (keeps within-trial correlation)."""
    ids = np.unique(groups[groups >= 0])
    members = {g: np.flatnonzero(groups == g) for g in ids}
    out = []
    for _ in range(nBoot):
        pick = rng.choice(ids, size=len(ids), replace=True)
        out.append(fn(np.concatenate([members[g] for g in pick])))
    return np.array(out)


def theta_transfer(d, trials, trialOf, rng):
    """theta = a0 + b u + g phi on force-on samples; design (G, G, -1)."""
    idx = np.concatenate([np.arange(tr['onStart'] + d['thetaLag'], tr['onEnd'] + 1)
                          for tr in trials])
    u = d['x'][idx] * d['c'][idx] + d['y'][idx] * d['s'][idx]
    A = np.column_stack((np.ones(len(idx)), u, d['ph'][idx]))
    yv = d['th'][idx]

    def fit(sel):
        return np.linalg.lstsq(A[sel], yv[sel], rcond=None)[0]

    co = fit(np.arange(len(idx)))
    boots = boot_by_trial(fit, trialOf[idx], N_BOOT, rng)
    design = np.array([THETA_GAIN, THETA_GAIN, -1.0])
    return dict(coef=co, err=boots.std(axis=0), design=design,
                rel=co / design - 1.0, r2=rsq(yv, A @ co),
                rms=float(np.std(yv - A @ co)),
                thetaMeas=yv, thetaDesign=A @ design)


def tangential(d, trialOf, rng):
    """First-order law r dpsi/dt = g theta (n.t) + F0.t on same-state sample pairs."""
    t, on = d['t'], d['on']
    i = np.arange(len(t) - 1)
    ok = (on[i] == on[i + 1]) & (np.diff(t) > 0) & (trialOf[i] >= 0)
    if d['fw'] < 2:
        ok &= np.r_[False, on[1:-1] == on[:-2]]
    i = i[ok]
    j = i + 1
    dpsi = wrap(d['psi'][j] - d['psi'][i])
    psiM = d['psi'][i] + 0.5 * dpsi
    tx, ty = -np.sin(psiM), np.cos(psiM)
    thM = 0.5 * (d['th'][i] + d['th'][j])
    A = np.column_stack((thM * (d['cEff'][i] * tx + d['sEff'][i] * ty), tx, ty))
    b = 0.5 * (d['r'][i] + d['r'][j]) * dpsi / (t[j] - t[i])

    def fit(sel):
        return np.linalg.lstsq(A[sel], b[sel], rcond=None)[0]

    co = fit(np.arange(len(b)))
    boots = boot_by_trial(fit, trialOf[i], N_BOOT, rng)
    onP = on[i]
    return dict(g=co[0], gErr=boots[:, 0].std(), F0=co[1:], F0err=boots[:, 1:].std(axis=0),
                r2on=rsq(b[onP], A[onP] @ co), A=A[onP], b=b[onP])


def predict_final(d, trials, g, F0, nSub=10):
    """Integrate the first-order law over each force-on window; final angle error (deg)."""
    t, th, r, psi = d['t'], d['th'], d['r'], d['psi']
    errs = []
    for tr in trials:
        i0 = tr['onStart'] + d['thetaLag']
        i1 = tr['onEnd']
        p = psi[i0]
        for i in range(i0, i1):
            h = (t[i + 1] - t[i]) / nSub
            for m in range(nSub):
                wgt = (m + 0.5) / nSub
                thW = (1 - wgt) * th[i] + wgt * th[i + 1]
                rW = (1 - wgt) * r[i] + wgt * r[i + 1]
                tx, ty = -math.sin(p), math.cos(p)
                p += h * (g * thW * (d['cEff'][i] * tx + d['sEff'][i] * ty)
                          + F0[0] * tx + F0[1] * ty) / rW
        errs.append(math.degrees(wrap(p - psi[i1])))
    return np.array(errs)


def probit_fit(m, o, rng, nBoot=N_BOOT):
    """ML fit of P(up) = Phi((m - mu)/sigma) with fractional outcomes o in {0, .5, 1}."""
    SIG_MIN, SIG_MAX, MU_MAX = 1e-3, 0.3, 0.15

    def grid_max(mm, oo, mus, lsig):
        MU, LS = np.meshgrid(mus, lsig, indexing='ij')
        p = normcdf((mm[None, None, :] - MU[..., None]) / np.exp(LS)[..., None])
        p = np.clip(p, 1e-12, 1 - 1e-12)
        ll = np.sum(oo * np.log(p) + (1 - oo) * np.log(1 - p), axis=-1)
        a, b = np.unravel_index(np.argmax(ll), ll.shape)
        return mus[a], lsig[b]

    def fit(mm, oo):
        # coarse grid, then a fine grid one coarse step around the best point
        mus = np.linspace(-MU_MAX, MU_MAX, 31)
        lsig = np.linspace(np.log(SIG_MIN), np.log(SIG_MAX), 31)
        m0, l0 = grid_max(mm, oo, mus, lsig)
        dm, dl = mus[1] - mus[0], lsig[1] - lsig[0]
        m1, l1 = grid_max(mm, oo, np.linspace(m0 - dm, m0 + dm, 21),
                          np.clip(np.linspace(l0 - dl, l0 + dl, 21),
                                  np.log(SIG_MIN), np.log(SIG_MAX)))
        return m1, float(np.exp(l1))

    mu, sg = fit(m, o)
    bs = np.array([fit(m[s], o[s]) for s in
                   (rng.integers(0, len(m), len(m)) for _ in range(min(nBoot, 200)))])
    return dict(mu=mu, sigma=sg, muCI=np.percentile(bs[:, 0], [16, 84]),
                sigmaCI=np.percentile(bs[:, 1], [16, 84]), atFloor=sg <= SIG_MIN * 1.01)


# =============================================================================
# GOAL 2: BORN RULE
# =============================================================================
def born_bins(T, probit):
    wCells = np.unique(np.round(T['w'], 4))
    h = float(np.median(np.diff(wCells))) if len(wCells) > 1 else 1.0
    U = np.linspace(-0.5 * h, 0.5 * h, 401)
    B = []
    for j in range(N_ANGLES):
        sel = T['bin'] == j
        n = int(sel.sum())
        if n == 0:
            continue
        o, P = T['outcome'][sel], T['P_born'][sel]
        # Each distinct phi value is one uniform draw w_c + U from its grid cell,
        # shared by every trial (repeat suite) that used it. Given U, trial i is up
        # with probability p_i(U) = Phi((P_i - w_c - U - mu)/sigma_m). Law of total
        # variance per cell: Var_U(sum p_i) -> hidden variable, E_U(sum p(1-p)) -> device.
        wr = np.round(T['w'][sel], 4)
        cells = np.unique(wr)
        varHV = varDev = upsIdeal = 0.0
        for wc in cells:
            upsIdeal += np.mean(np.sum((wc + U[:, None]) < P[wr == wc][None, :], axis=1))
            pU = normcdf((P[wr == wc][None, :] - wc - U[:, None] - probit['mu'])
                         / probit['sigma'])
            varHV += np.var(pU.sum(axis=1))
            varDev += np.mean(np.sum(pU * (1 - pU), axis=1))
        sHV = math.sqrt(varHV) / n
        sDev = math.sqrt(varDev) / n
        nAmb = int(T['ambiguous'][sel].sum())
        sAmb = math.sqrt(0.25 * nAmb) / n
        B.append(dict(bin=j, D=float(np.mean(T['D'][sel])), n=n, nPhi=len(cells),
                      nAmb=nAmb, p=float(o.mean()), P=float(P.mean()),
                      E=float(T['ideal_up'][sel].mean()),
                      gridBias=upsIdeal / n - float(P.mean()),
                      sHV=sHV, sDev=sDev, sAmb=sAmb,
                      sig=max(math.sqrt(sHV ** 2 + sDev ** 2 + sAmb ** 2), SIGMA_FLOOR)))
    B = {k: np.array([b[k] for b in B]) for k in B[0]}
    chi2 = float(np.sum(((B['p'] - B['P']) / B['sig']) ** 2))
    dof = len(B['p'])
    return B, chi2, dof, chi2_sf(chi2, dof), h


# =============================================================================
# FIGURES
# =============================================================================
def fig_trial_cycle(d, trials, path):
    a, b = trials[TRIALCYCLE[0]]['start'], trials[TRIALCYCLE[1]]['end']
    sl = slice(a, b + 1)
    t, on = d['t'][sl], d['on'][sl]
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(6.4, 4.8))
    ax[0].plot(t, d['c'][sl] * on, label='c')
    ax[0].plot(t, d['s'][sl] * on, label='s')
    ax[0].plot(t, d['ph'][sl], label=r'$\phi$')
    ax[1].plot(t, d['x'][sl], '--', label='x')
    ax[1].plot(t, d['y'][sl], '--', label='y')
    ax[1].plot(t, d['th'][sl], '--', label=r'$\theta$')
    for axx in ax:
        axx.fill_between(t, 0, 1, where=on, color='0.9', transform=axx.get_xaxis_transform(),
                         zorder=0, lw=0)
        axx.legend(loc='upper right', fontsize=8)
    ax[0].set_ylabel('Control signals')
    ax[1].set_ylabel('Measurements')
    ax[1].set_xlabel('t (s)')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_empirical(B, chi2, dof, pval, path):
    fig, (ax, axr) = plt.subplots(2, 1, sharex=True, figsize=(6.4, 5.4),
                                  gridspec_kw={'height_ratios': [3, 1]})
    Dg = np.linspace(0, ANGLE_MAX, 400)
    ax.plot(Dg, np.cos(0.5 * Dg) ** 2, color='tab:blue', label='Born rule', zorder=1)
    ax.plot(B['D'], B['E'], 's', mfc='none', mec='0.45', ms=8,
            label=r'ideal model, same $\phi$ values', zorder=2)
    ax.errorbar(B['D'], B['p'], yerr=B['sig'], fmt='.', color='red', ecolor='k',
                ms=14, capsize=4, label='empirical', zorder=3)
    ax.set_ylabel('Up probability')
    ax.legend(fontsize=8, loc='lower left')
    ax.set_title(rf'$\chi^2/\nu = {chi2:.1f}/{dof}$,  $p = {pval:.2f}$', fontsize=10)
    pull = (B['p'] - B['P']) / B['sig']
    axr.axhspan(-1, 1, color='0.9', lw=0)
    axr.axhline(0, color='k', lw=0.8)
    axr.plot(B['D'], pull, 'o', color='red', ms=5)
    axr.set_ylim(-3.5, 3.5)
    axr.set_ylabel(r'pull ($\sigma$)')
    axr.set_xlabel('Relative angle (rad)')
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fig_ode(TH, TG, T, probit, errDesign, path):
    fig, axs = plt.subplots(2, 2, figsize=(9, 8))
    # (a) theta transfer
    ax = axs[0, 0]
    ax.plot(TH['thetaDesign'], TH['thetaMeas'], '.', ms=2, alpha=0.3)
    lim = [min(TH['thetaDesign'].min(), TH['thetaMeas'].min()),
           max(TH['thetaDesign'].max(), TH['thetaMeas'].max())]
    ax.plot(lim, lim, 'k--', lw=1)
    ax.set_xlabel(r'designed $\theta = G(1+u) - \phi$  (/5 V)')
    ax.set_ylabel(r'logged $\theta$  (/5 V)')
    ax.set_title(f'(a) theta transfer, rms = {TH["rms"]:.4f}', fontsize=10)
    # (b) tangential dynamics
    ax = axs[0, 1]
    pred = TG['A'] @ np.r_[TG['g'], TG['F0']]
    ax.plot(pred, TG['b'], '.', ms=2, alpha=0.3)
    lim = np.percentile(np.abs(TG['b']), 99.5)
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel(r'model $g\,\theta\,(\hat n\cdot\hat t) + F_0\cdot\hat t$  (/s)')
    ax.set_ylabel(r'measured $r\,\dot\psi$  (/s)')
    ax.set_title(rf'(b) first-order law, g = {TG["g"]:.1f} $\pm$ {TG["gErr"]:.1f} /s, '
                 rf'$R^2$ = {TG["r2on"]:.3f}', fontsize=10)
    # (c) outcome vs margin with threshold fit
    ax = axs[1, 0]
    col = np.where(T['model_ok'], 'tab:green', 'tab:red')
    jit = np.random.default_rng(1).uniform(-0.03, 0.03, len(T['margin']))
    ax.scatter(T['margin'], T['outcome'] + jit, c=col, s=10)
    mg = np.linspace(T['margin'].min(), T['margin'].max(), 400)
    ax.plot(mg, normcdf((mg - probit['mu']) / probit['sigma']), 'k-', lw=1)
    ax.axvline(0, color='0.5', lw=0.8, ls=':')
    ins = ax.inset_axes([0.55, 0.12, 0.4, 0.45])
    zoom = np.abs(T['margin']) < 0.06
    ins.scatter(T['margin'][zoom], T['outcome'][zoom] + jit[zoom], c=col[zoom], s=10)
    mz = np.linspace(-0.06, 0.06, 200)
    ins.plot(mz, normcdf((mz - probit['mu']) / probit['sigma']), 'k-', lw=1)
    ins.axvline(0, color='0.5', lw=0.8, ls=':')
    ins.set_xlim(-0.06, 0.06)
    ins.tick_params(labelsize=7)
    ax.set_xlabel(r'margin $m = \cos^2(D/2) - w$  (ideal model: up iff $m>0$)')
    ax.set_ylabel('outcome (1 = up)')
    ax.set_title(rf'(c) threshold: $\mu$ = {probit["mu"]:+.3f}, '
                 rf'$\sigma_m$ = {probit["sigma"]:.3f}', fontsize=10)
    # (d) forward prediction with design g
    ax = axs[1, 1]
    lim = max(2.0, float(np.ceil(np.percentile(np.abs(errDesign), 99))))
    ax.hist(np.clip(errDesign, -lim, lim), bins=np.linspace(-lim, lim, 41),
            color='tab:blue', alpha=0.8)
    ax.set_xlabel(r'final angle, predicted $-$ measured (deg)')
    ax.set_ylabel('trials')
    ax.set_title(rf'(d) forward prediction, design g: median |err| = '
                 rf'{np.median(np.abs(errDesign)):.1f}$^\circ$', fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


# =============================================================================
# GOAL 3: ANIMATIONS
# =============================================================================
def save_anim(anim, path, fps):
    if path.endswith('.mp4'):
        anim.save(path, writer=ani.FFMpegWriter(fps=fps, bitrate=2400), dpi=150)
    else:
        anim.save(path, writer=ani.PillowWriter(fps=fps), dpi=100)


def anim_trial(d, trials, T, k, path, fps):
    tr = trials[k]
    a = max(0, tr['pre'] - 8)
    b = tr['end']
    frames = np.arange(a, b + 1)
    t0 = d['t'][tr['onStart']]
    R = float(np.median(d['r']))
    axis = math.atan2(d['s'][tr['onStart']], d['c'][tr['onStart']])
    ca, sa = math.cos(axis), math.sin(axis)

    fig = plt.figure(figsize=(9, 4.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25])
    axP = fig.add_subplot(gs[:, 0], aspect='equal')
    axT = fig.add_subplot(gs[0, 1])
    axA = fig.add_subplot(gs[1, 1], sharex=axT)

    # plane, natural units (logged / median radius)
    axP.add_patch(Circle((0, 0), 1, fill=False, ec='0.6'))
    axP.plot([0, ca], [0, sa], color='orange', lw=2)
    axP.plot([0, -ca], [0, -sa], color='orange', lw=2, ls='--')
    axP.plot(ca, sa, 'o', color='orange', mfc='white', ms=9)
    axP.text(1.12 * ca, 1.12 * sa, r'$\uparrow$', ha='center', va='center', fontsize=14)
    axP.text(-1.12 * ca, -1.12 * sa, r'$\downarrow$', ha='center', va='center', fontsize=14)
    ip = tr['pre']
    axP.plot(d['x'][ip] / R, d['y'][ip] / R, 'x', color='0.4', ms=8)
    axP.set_xlim(-1.3, 1.3); axP.set_ylim(-1.3, 1.3)
    axP.set_xticks([]); axP.set_yticks([])
    trail, = axP.plot([], [], '-', color='tab:blue', alpha=0.35, lw=2)
    dot, = axP.plot([], [], 'o', color='tab:blue', ms=10)
    lamp = axP.add_patch(Circle((1.12, 1.12), 0.07, color='red'))
    axP.text(1.12, 0.98, 'force on', ha='center', va='top', fontsize=7, color='0.4')

    # time series
    tt = d['t'][frames] - t0
    onF = d['on'][frames]
    for axx in (axT, axA):
        axx.fill_between(tt, 0, 1, where=onF, color='0.92', lw=0,
                         transform=axx.get_xaxis_transform())
    axT.plot(tt, d['th'][frames], color='tab:green', label=r'$\theta$')
    axT.plot(tt, d['ph'][frames], color='tab:purple', lw=1, label=r'$\phi$')
    axT.axhline(0, color='k', lw=0.6)
    axT.legend(fontsize=8, loc='upper right')
    axT.set_ylabel('/5 V')
    rel = np.degrees(np.unwrap(wrap(d['psi'][frames] - axis)))
    axA.plot(tt, rel, color='tab:blue')
    for yl in (0, 180, -180):
        axA.axhline(yl, color='orange', lw=0.8, ls='-' if yl == 0 else '--')
    axA.set_ylabel(r'$\psi - \psi_{axis}$ (deg)')
    axA.set_xlabel('time since force on (s)')
    curT = axT.axvline(tt[0], color='k', lw=0.8)
    curA = axA.axvline(tt[0], color='k', lw=0.8)

    lab = {1.0: 'up', 0.0: 'down', 0.5: 'ambiguous'}
    ok = 'OK' if (T['model_ok'][k] and T['sign_ok'][k]) else 'MISMATCH'
    fig.suptitle(f'Trial {k}:  D = {T["D_deg"][k]:.0f} deg,  w(phi) = {T["w"][k]:.2f},  '
                 f'P_Born = {T["P_born"][k]:.2f},  theta0 = {T["theta0"][k]:+.3f}\n'
                 f'outcome {lab[T["outcome"][k]]},  ideal model {lab[T["ideal_up"][k]]}'
                 f'  [{ok}]', fontsize=10)
    fig.tight_layout()

    def step(n):
        i = frames[n]
        lo = max(0, n - 10)
        trail.set_data(d['x'][frames[lo:n + 1]] / R, d['y'][frames[lo:n + 1]] / R)
        dot.set_data([d['x'][i] / R], [d['y'][i] / R])
        lamp.set_visible(bool(d['on'][i]))
        curT.set_xdata([tt[n], tt[n]])
        curA.set_xdata([tt[n], tt[n]])
        return trail, dot, lamp, curT, curA

    save_anim(ani.FuncAnimation(fig, step, frames=len(frames)), path, fps)
    plt.close(fig)


def anim_overview(d, trials, T, path, fps, nPre=4):
    """All trials at once in the axis frame: up pole at top, down at bottom."""
    nOn = min(tr['onEnd'] - tr['onStart'] + 1 for tr in trials)
    lags = np.arange(-nPre, nOn)
    rel = np.empty((len(trials), len(lags)))
    for k, tr in enumerate(trials):
        axis = math.atan2(d['s'][tr['onStart']], d['c'][tr['onStart']])
        for n, L in enumerate(lags):
            i = tr['pre'] + 1 + L if L < 0 else tr['onStart'] + L
            rel[k, n] = wrap(d['psi'][i] - axis)
    dt = float(np.median(np.diff(d['t'])))
    bins = T['bin']
    Dc = np.arange(N_ANGLES) * ANGLE_MAX / (N_ANGLES - 1)
    cmap = plt.get_cmap('viridis')

    fig, (axC, axB) = plt.subplots(1, 2, figsize=(9, 4.4),
                                   gridspec_kw={'width_ratios': [1, 1.2]})
    axC.set_aspect('equal')
    axC.add_patch(Circle((0, 0), 1, fill=False, ec='0.6'))
    axC.plot([0, 0], [0, 1], color='orange', lw=2)
    axC.plot([0, 0], [0, -1], color='orange', lw=2, ls='--')
    axC.text(0, 1.15, r'$\uparrow$', ha='center', fontsize=14)
    axC.text(0, -1.25, r'$\downarrow$', ha='center', fontsize=14)
    axC.set_xlim(-1.3, 1.3); axC.set_ylim(-1.35, 1.3)
    axC.axis('off')
    rng = np.random.default_rng(2)
    rad = 1.0 + rng.uniform(-0.05, 0.05, len(trials))   # jitter so stacked dots show
    sc = axC.scatter(np.zeros(len(trials)), np.zeros(len(trials)), s=14,
                     c=cmap(bins / (N_ANGLES - 1)), edgecolors='none')
    clock = axC.text(-1.3, -1.35, '', fontsize=9)

    Dg = np.linspace(0, ANGLE_MAX, 200)
    axB.plot(np.degrees(Dg), np.cos(0.5 * Dg) ** 2, color='tab:blue', label='Born rule')
    bars = axB.bar(np.degrees(Dc), np.zeros(N_ANGLES), width=7,
                   color=cmap(np.arange(N_ANGLES) / (N_ANGLES - 1)), alpha=0.8,
                   label='fraction nearer up pole')
    axB.set_ylim(0, 1.05)
    axB.set_xlabel('relative angle D (deg)')
    axB.set_ylabel('fraction')
    axB.legend(fontsize=8, loc='lower left')
    fig.suptitle(f'{len(trials)} measurements in the measurement-axis frame', fontsize=11)
    fig.tight_layout()

    def step(n):
        ang = np.pi / 2 + rel[:, n]
        sc.set_offsets(np.column_stack((rad * np.cos(ang), rad * np.sin(ang))))
        upNow = np.cos(rel[:, n]) > 0
        for j, bar in enumerate(bars):
            sel = bins == j
            bar.set_height(upNow[sel].mean() if sel.any() else 0)
        state = 'force off' if lags[n] < 0 else 'force on'
        clock.set_text(f't = {lags[n] * dt:+.2f} s ({state})')
        return (sc, clock, *bars)

    hold = [len(lags) - 1] * int(1.5 * fps)          # linger on the final frame
    save_anim(ani.FuncAnimation(fig, step, frames=list(range(len(lags))) + hold), path, fps)
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================
def main(logfile=LOGFILE, outdir=OUTDIR, animate=ANIMATE, fmt=ANIM_FORMAT):
    rng = np.random.default_rng(SEED)
    here = os.path.dirname(os.path.abspath(__file__))
    logPath = logfile if os.path.isabs(logfile) else os.path.join(here, logfile)
    stem = os.path.splitext(os.path.basename(logPath))[0]
    outdir = outdir or os.path.join(here, 'analysis_' + stem)
    os.makedirs(outdir, exist_ok=True)
    lines = []

    def out(s=''):
        print(s)
        lines.append(s)

    d = load_log(logPath)
    trials, trialOf = segment(d)
    T = trial_table(d, trials)
    out(f'Log: {os.path.basename(logPath)}  (fw v{d["fw"]}, {len(d["t"])} samples, '
        f'{len(trials)} trials, median dt = {np.median(np.diff(d["t"])) * 1e3:.0f} ms)')

    # ---------------- GOAL 1 ----------------
    TH = theta_transfer(d, trials, trialOf, rng)
    Robs = float(np.median(d['r'][trialOf >= 0]))
    Riqr = np.percentile(d['r'][trialOf >= 0], [25, 75])
    TG = tangential(d, trialOf, rng)
    errFit = predict_final(d, trials, TG['g'], TG['F0'])
    errDesign = predict_final(d, trials, G_DESIGN, (0.0, 0.0))
    probit = probit_fit(T['margin'], T['outcome'], rng)
    signRate = float(T['sign_ok'].mean())
    agreeRate = float(T['model_ok'].mean())

    out('\n=== GOAL 1: does the circuit implement the ODE model? ===')
    out('1a theta transfer  theta = a0 + b u + g phi   (fit +/- bootstrap, design, deviation)')
    for name, c_, e_, ds, rl in zip(('a0', 'b', 'g'), TH['coef'], TH['err'],
                                    TH['design'], TH['rel']):
        out(f'     {name:2s} = {c_:+.4f} +/- {e_:.4f}   design {ds:+.4f}   {100 * rl:+.1f}%')
    out(f'     R^2 = {TH["r2"]:.4f}, rms residual = {TH["rms"]:.4f} (/5 V)')
    out(f'1b radius  R_obs = {Robs:.4f} (IQR {Riqr[0]:.3f}-{Riqr[1]:.3f})   '
        f'design sqrt(14/19) = {R_EQ:.4f}   {100 * (Robs / R_EQ - 1):+.1f}%')
    out(f'1c dynamics  r dpsi/dt = g theta (n.t) + F0.t')
    out(f'     g = {TG["g"]:.1f} +/- {TG["gErr"]:.1f} /s   design {G_DESIGN:.0f} /s   '
        f'{100 * (TG["g"] / G_DESIGN - 1):+.1f}%   (R^2 force-on = {TG["r2on"]:.3f})')
    out(f'     drift F0 = ({TG["F0"][0]:+.4f} +/- {TG["F0err"][0]:.4f}, '
        f'{TG["F0"][1]:+.4f} +/- {TG["F0err"][1]:.4f}) /s')
    out(f'     forward prediction of final angle: median |err| = '
        f'{np.median(np.abs(errDesign)):.1f} deg (design g, no drift), '
        f'{np.median(np.abs(errFit)):.1f} deg (fitted g, F0)')
    out(f'1d sign rule  pole = sign(theta0)*(c,s): {int(T["sign_ok"].sum())}/{len(T["trial"])}'
        f' = {100 * signRate:.1f}%')
    bad = np.flatnonzero(~T['sign_ok'])
    if len(bad):
        out(f'     exceptions (trial: theta0): '
            + ', '.join(f'{k}: {T["theta0"][k]:+.3f}' for k in bad))
    out(f'1e end to end  outcome == ideal model(D, phi): {int(T["model_ok"].sum())}/'
        f'{len(T["trial"])} = {100 * agreeRate:.1f}%   ({int(T["ambiguous"].sum())} ambiguous)')
    sgTxt = (f'< {probit["sigmaCI"][1]:.3f} (no outcomes disagree near threshold)'
             if probit['atFloor'] else
             f'{probit["sigma"]:.3f} (68% CI {probit["sigmaCI"][0]:.3f}-{probit["sigmaCI"][1]:.3f})')
    out(f'     threshold bias mu = {probit["mu"]:+.3f} '
        f'(68% CI {probit["muCI"][0]:+.3f} to {probit["muCI"][1]:+.3f}),  '
        f'jitter sigma_m = {sgTxt}   [probability units]')
    mis = np.flatnonzero(~T['model_ok'])
    if len(mis):
        out('     mismatches (trial: margin m): '
            + ', '.join(f'{k}: {T["margin"][k]:+.3f}' for k in mis))

    # ---------------- GOAL 2 ----------------
    B, chi2, dof, pval, h = born_bins(T, probit)
    out('\n=== GOAL 2: do the empirical probabilities obey the Born rule? ===')
    out(f'phi grid: {int(np.max(B["nPhi"]))} distinct values, spacing h = {h:.3f} in w '
        f'-> hidden-variable resolution <= 1/(2 N_phi) = {0.5 / np.max(B["nPhi"]):.3f}')
    out(' bin  D(deg)   n  amb   p_hat   Born   ideal  gridBias  sig_HV  sig_dev sig_amb  sigma   pull')
    for q in range(len(B['p'])):
        out(f' {B["bin"][q]:3d}  {round(np.degrees(B["D"][q]), 1) + 0.0:6.1f}  {B["n"][q]:3d}  {B["nAmb"][q]:3d}'
            f'   {B["p"][q]:.3f}  {B["P"][q]:.3f}  {B["E"][q]:.3f}  {B["gridBias"][q]:+.3f}'
            f'   {B["sHV"][q]:.3f}   {B["sDev"][q]:.3f}   {B["sAmb"][q]:.3f}'
            f'   {B["sig"][q]:.3f}  {(B["p"][q] - B["P"][q]) / B["sig"][q]:+.2f}')
    out(f'chi^2 = {chi2:.2f} for {dof} dof (chi^2/nu = {chi2 / dof:.2f}), p = {pval:.3f}')
    out(f'RMS(p_hat - Born) = {np.sqrt(np.mean((B["p"] - B["P"]) ** 2)):.3f};  '
        f'RMS(ideal model on this grid - Born) = {np.sqrt(np.mean((B["E"] - B["P"]) ** 2)):.3f}')
    out('ideal    = what a perfect circuit gives with exactly these phi values')
    out('gridBias = expected p_hat - Born for a perfect circuit under this phi grid '
        '(nonzero for an endpoint-inclusive grid; zero for a midpoint grid)')

    # ---------------- VERDICT ----------------
    def grade(key, val):
        lo, hi = TOL[key]['max']
        if not lo <= val <= hi:
            return 'FAIL'
        lo, hi = TOL[key]['typ']
        return 'PASS' if lo <= val <= hi else 'WARN'

    def band(key, pct=True):
        f = (lambda v: f'{100 * v:+.1f}%') if pct else (lambda v: f'{v:+.3f}')
        (a1, b1), (a2, b2) = TOL[key]['typ'], TOL[key]['max']
        return f'typ [{f(a1)}, {f(b1)}], max [{f(a2)}, {f(b2)}]'

    checks = []
    for key, val in zip(('a0', 'b', 'gphi'), TH['rel']):
        checks.append((grade(key, val), f'1a theta {key}', f'{100 * val:+.1f}%', band(key)))
    relR, relG = Robs / R_EQ - 1, TG['g'] / G_DESIGN - 1
    checks.append((grade('R', relR), '1b equilibrium radius', f'{100 * relR:+.1f}%', band('R')))
    checks.append((grade('g', relG), '1c measurement gain g', f'{100 * relG:+.1f}%', band('g')))
    checks.append((grade('mu', probit['mu']), '1e threshold bias mu', f'{probit["mu"]:+.3f}',
                   band('mu', pct=False)))
    checks.append(('PASS' if pval >= P_VALUE_MIN else 'FAIL', '2  Born rule chi^2',
                   f'p = {pval:.3f}', f'need p >= {P_VALUE_MIN}'))
    out('\n=== SUMMARY ===  (WARN: outside the typical-offset band; FAIL: outside the max band)')
    for status, name, val, rng_ in checks:
        out(f'  [{status}] {name:24s} {val:>8s}   {rng_}')
    out(f'  [info] 1d sign rule {100 * signRate:.1f}%, 1e ideal-model agreement '
        f'{100 * agreeRate:.1f}% (outcome-level agreement is tested by the chi^2)')

    # ---------------- FILES ----------------
    with open(os.path.join(outdir, 'summary.txt'), 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')
    keys = ['trial', 't0', 'D_deg', 'bin', 'phi', 'w', 'P_born', 'margin', 'theta0',
            'final_deg', 'outcome', 'ideal_up', 'sign_ok', 'model_ok']
    with open(os.path.join(outdir, 'trials.csv'), 'w', newline='') as f:
        wr = csv.writer(f)
        wr.writerow(keys)
        for n in range(len(T['trial'])):
            wr.writerow([f'{T[k][n]:.4f}' if isinstance(T[k][n], float) else int(T[k][n])
                         for k in keys])
    fig_trial_cycle(d, trials, os.path.join(outdir, 'trialCycle.pdf'))
    fig_empirical(B, chi2, dof, pval, os.path.join(outdir, 'empiricalProb.pdf'))
    fig_ode(TH, TG, T, probit, errDesign, os.path.join(outdir, 'odeChecks.pdf'))
    print(f'\nFigures, summary.txt and trials.csv written to {outdir}')

    # ---------------- ANIMATIONS ----------------
    if fmt == 'mp4' and not ani.writers.is_available('ffmpeg'):
        print('ffmpeg not found; writing GIFs instead')
        fmt = 'gif'
    fps = max(1, int(round(ANIM_SPEED / np.median(np.diff(d['t'])))))
    animDir = os.path.join(outdir, 'animations')
    os.makedirs(animDir, exist_ok=True)
    if animate == 'all':
        which = list(range(len(trials)))
    elif animate == 'failures':
        which = list(np.flatnonzero(~(T['model_ok'] & T['sign_ok'])))
    elif animate in (None, 'none'):
        which = []
    else:
        which = [int(k) for k in animate]
    for k in which:
        path = os.path.join(animDir, f'trial_{k:03d}.{fmt}')
        if os.path.exists(path) and not OVERWRITE:
            continue
        print(f'  animating trial {k}')
        anim_trial(d, trials, T, k, path, fps)
    if OVERVIEW:
        path = os.path.join(animDir, f'overview.{fmt}')
        if OVERWRITE or not os.path.exists(path):
            print('  animating overview')
            anim_overview(d, trials, T, path, max(1, fps // 2))
    return dict(d=d, trials=trials, T=T, TH=TH, TG=TG, probit=probit, B=B,
                chi2=chi2, pval=pval, checks=checks)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    ap.add_argument('log', nargs='?', default=LOGFILE)
    ap.add_argument('--out', default=OUTDIR)
    ap.add_argument('--anim', default=ANIMATE,
                    help="none | failures | all | comma-separated trial numbers")
    ap.add_argument('--fmt', default=ANIM_FORMAT, choices=('gif', 'mp4'))
    args, _ = ap.parse_known_args()      # tolerate Spyder/IPython extra args
    anim = args.anim
    if isinstance(anim, str) and anim not in ('none', 'failures', 'all'):
        anim = [int(v) for v in anim.split(',')]
    res = main(args.log, args.out, anim, args.fmt)
