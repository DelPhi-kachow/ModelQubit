# -*- coding: utf-8 -*-
"""
toleranceMC.py -- propagate component tolerances through the model-qubit circuit
to the quantities mqAnalysis.py checks, in the same (logged, /5 V) units.

Circuits (DynamicalCollapsePaper/model_qubit_circuit_diagrams):
  theta.tex    two AD633 (x c, y s) -> 22k each into the + node, which also has
               22k+22k from +5 V and 10k to ground; phi via 10k into -, 10k feedback
  thetacs.tex  AD633: theta * c (theta on X, c on Y)
  xy.tex       velocity stage: 1k in (theta c), 1k in (radial), 1k || 10 nF feedback;
               position integrator 1 MOhm, 10 nF; AD633 squarer x*x; AD633 x*(1-r^2)
  rsquare.tex  x^2, y^2 via 10k into -, 100k feedback; + from 47k/10k divider off +5 V
  ardSK/ard1   c/s drive: DC-unity filter, then gain (1 + Rf/Rg) with Rg to +5 V
  phi drive    22k/22k divider, unity follower (ard2_frag.tex not yet updated)
  ard3.tex     theta readout: + node gets theta via 10k, +5 V via 10k+10k, 10k to
               ground; gain 1 + 10k/10k  ->  V_pin = 2 V + 0.8 V_theta
  ard4.tex     x/y readout:  + node gets x via 10k, +5 V via 10k, 10k to ground;
               gain 1 + (10k||10k)/10k  ->  V_pin = 2.5 V + 0.5 V_x
The Arduino ADC reference and all +5 V references are taken as one rail V_S.

Tolerances: resistors uniform +/-1%, capacitors uniform +/-10%, AD633 scale
factor uniform +/-0.25%. AD633 offsets (datasheet): output +/-5 mV typ / +/-50 mV
max, X/Y input +/-5 mV typ / +/-30 mV max, run as two scenarios.
"""
import numpy as np

N = 40000
rng = np.random.default_rng(0)
VS_NOM = 5.0


def R(nom, tol=0.01):
    return nom * (1 + rng.uniform(-tol, tol, N))


def ad633(offOut, offIn, scaleTol=0.0025):
    k = 1 + rng.uniform(-scaleTol, scaleTol, N)
    oo = rng.uniform(-offOut, offOut, N)
    ox = rng.uniform(-offIn, offIn, N)
    oy = rng.uniform(-offIn, offIn, N)
    return lambda X, Y: k * (X + ox) * (Y + oy) / 10.0 + oo, k, ox


def par(a, b):
    return a * b / (a + b)


def simulate(offOut, offIn, VS=VS_NOM, resTol=0.01, capTol=0.10, scaleTol=0.0025):
    Rn = lambda nom: R(nom, resTol)
    AD = lambda: ad633(offOut, offIn, scaleTol)
    # ---- readout (ard3, ard4) and drive (ardSK, phi) ----
    th_a, th_b, th_c, th_f, th_g = Rn(10e3), Rn(10e3) + Rn(10e3), Rn(10e3), Rn(10e3), Rn(10e3)
    xr = [dict(a=Rn(10e3), b=Rn(10e3), c=Rn(10e3), f=par(Rn(10e3), Rn(10e3)), g=Rn(10e3))
          for _ in range(2)]
    cs = [dict(f=Rn(10e3), g=Rn(10e3)) for _ in range(2)]
    ph1, ph2 = Rn(22e3), Rn(22e3)

    def theta_log(Vth):
        vp = (Vth / th_a + VS / th_b) / (1 / th_a + 1 / th_b + 1 / th_c)
        return 1.25 * (1 + th_f / th_g) * vp / VS - 0.5

    def q_log(Vq, ch):
        r = xr[ch]
        vp = (Vq / r['a'] + VS / r['b']) / (1 / r['a'] + 1 / r['b'] + 1 / r['c'])
        return 2 * (1 + r['f'] / r['g']) * vp / VS - 1

    def drive(cf, ch):
        d = cs[ch]
        vin = (cf + 1) / 2 * VS
        return vin * (1 + d['f'] / d['g']) - VS * d['f'] / d['g']

    def phi_V(pf):
        return 2 * pf * VS * ph2 / (ph1 + ph2)

    # ---- theta circuit ----
    mXC, _, _ = AD()
    mYS, _, _ = AD()
    w1, w2, t1, t2, tg, ti, tf = (Rn(22e3), Rn(22e3), Rn(22e3), Rn(22e3), Rn(10e3),
                                  Rn(10e3), Rn(10e3))

    def theta_V(Vx, Vy, Vc, Vs_, Vphi):
        W1, W2 = mXC(Vx, Vc), mYS(Vy, Vs_)
        G = 1 / w1 + 1 / w2 + 1 / (t1 + t2) + 1 / tg
        vp = (W1 / w1 + W2 / w2 + VS / (t1 + t2)) / G
        return (1 + tf / ti) * vp - (tf / ti) * Vphi

    # ---- radius (rsquare + squarers + radial multiplier), per channel ----
    rf, r1, r2, rt, rgd = Rn(100e3), Rn(10e3), Rn(10e3), Rn(47e3), Rn(10e3)
    sq = [AD() for _ in range(2)]
    radM = [AD() for _ in range(2)]
    Vplus = VS * rgd / (rt + rgd)
    radii = []
    for ch in range(2):
        other = sq[1 - ch][0](0.0, 0.0)          # other squarer with its input at 0
        o_rad = radM[ch][2]                      # (1-r^2) enters the radial AD633 on X
        # equilibrium: V_{1-r^2} = -o_rad (radial force x*(V+o)/10 = 0)
        target = (1 + rf / par(r1, r2)) * Vplus - (rf / r2 if ch == 0 else rf / r1) * other + o_rad
        # solve (rf/r_ch) * sq(V, V) = target for V by bisection on [3, 5.5] V
        lo, hi = np.full(N, 3.0), np.full(N, 5.5)
        gain = rf / r1 if ch == 0 else rf / r2
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            f = gain * sq[ch][0](mid, mid) - target
            lo, hi = np.where(f > 0, lo, mid), np.where(f > 0, mid, hi)
        V = 0.5 * (lo + hi)
        radii.append((q_log(V, ch) - q_log(-V, ch)) / 2)   # symmetric readout (offset cancels)
    Rlog = 0.5 * (radii[0] + radii[1])
    Rvolt = VS_NOM * 0.5 * (radii[0] + radii[1])  # only for the geometry below

    # ---- theta transfer fit, as mqAnalysis does it ----
    # Sample states on the circle and axes; build logged (u, phi) -> logged theta
    ang = np.linspace(0, 2 * np.pi, 9)[:-1]
    D = np.linspace(0, np.pi, 7)
    phis = np.linspace(0.03, 0.40, 5)
    rows, ys = [], []
    for a in ang:
        for dd in D:
            for pf in phis:
                cf, sf = np.cos(a + dd), np.sin(a + dd)
                Vx, Vy = 4.29 * np.cos(a), 4.29 * np.sin(a)
                Vth = theta_V(Vx, Vy, drive(cf, 0), drive(sf, 1), phi_V(pf))
                u = q_log(Vx, 0) * cf + q_log(Vy, 1) * sf
                rows.append((u, np.full(N, pf)))
                ys.append(theta_log(Vth))
    U = np.array([r[0] for r in rows])          # (M, N)
    P = np.array([r[1] for r in rows])
    Y = np.array(ys)
    coef = np.empty((N, 3))
    for n in range(N):                          # small lstsq per sample
        A = np.column_stack((np.ones(len(U)), U[:, n], P[:, n]))
        coef[n] = np.linalg.lstsq(A, Y[:, n], rcond=None)[0]

    # ---- measurement gain g (first-order law) in logged units ----
    mTC, kTC, oxTC = AD()
    rfv, rin, rpos = Rn(1e3), Rn(1e3), Rn(1e6)
    cpos = 10e-9 * (1 + rng.uniform(-capTol, capTol, N))
    # d x_log / dt per unit (theta_log * c_float), small-signal around theta = 0.2
    th0, eps = 0.2, 1e-3
    Vth = lambda tl: (tl + 0.5) / 1.25 * VS  # nominal inverse is fine for a derivative
    dVth = (Vth(th0 + eps) - Vth(th0)) / (theta_log(Vth(th0 + eps)) - theta_log(Vth(th0)))
    dVc = drive(1.0, 0) - drive(0.0, 0)
    dxlog_dVx = (q_log(0.1, 0) - q_log(0.0, 0)) / 0.1
    g = dxlog_dVx * (rfv / rin) * kTC / (10 * rpos * cpos) * dVth * dVc * eps / eps
    # ^ dVth already converts one unit of theta_log to volts

    # ---- Born threshold bias mu(D) = cos^2(D/2) - P_circuit(D) ----
    PHI_MIN = 10 / 47 * (1 - np.sqrt(14 / 19))
    PHI_MAX = 10 / 47 * (1 + np.sqrt(14 / 19))
    Dg = np.linspace(0, np.pi / 2, 10)
    mu = np.empty((N, len(Dg)))
    for j, dd in enumerate(Dg):
        cf, sf = np.cos(dd), np.sin(dd)
        Vx, Vy = Rvolt, 0.0 * Rvolt
        f = lambda pf: theta_V(Vx, Vy, drive(cf, 0), drive(sf, 1), phi_V(pf)) + oxTC
        f0, f1 = f(0.0), f(0.5)
        pstar = -f0 / ((f1 - f0) / 0.5)          # theta_eff linear in phi
        w = (pstar - PHI_MIN) / (PHI_MAX - PHI_MIN)
        mu[:, j] = np.cos(dd / 2) ** 2 - np.clip(w, 0, 1)
    return dict(a0=coef[:, 0], b=coef[:, 1], gphi=coef[:, 2], R=Rlog, g=g, mu=mu)


DESIGN = dict(a0=10 / 47, b=10 / 47, gphi=-1.0, R=np.sqrt(14 / 19), g=50.0)
MEAS = dict(a0=0.2119, b=0.2096, gphi=-1.0022, R=0.8655, g=45.5, mu=-0.005)

if __name__ == '__main__':
    nom = simulate(0, 0, resTol=0, capTol=0, scaleTol=0)
    print('nominal check:', {k: float(np.mean(v)) for k, v in nom.items() if k != 'mu'},
          'mu', float(np.abs(nom['mu']).max()))
    print('\nRelative deviation from design (%): central 95% interval, and measured')
    res = {}
    for label, oo, oi in (('typ offsets', 5e-3, 5e-3), ('max offsets', 50e-3, 30e-3)):
        res[label] = simulate(oo, oi)
    for k in ('a0', 'b', 'gphi', 'R', 'g'):
        line = f'  {k:5s}'
        for label, S in res.items():
            rel = 100 * (S[k] / DESIGN[k] - 1)
            lo, hi = np.percentile(rel, [2.5, 97.5])
            line += f'   {label}: [{lo:+6.2f}, {hi:+6.2f}]'
        line += f'   measured {100 * (MEAS[k] / DESIGN[k] - 1):+6.2f}'
        print(line)
    for label, S in res.items():
        mbar = S['mu'].mean(axis=1)
        mmax = S['mu'][np.arange(N), np.argmax(np.abs(S['mu']), axis=1)]
        print(f'  mu (prob. units), {label}: mean over D 95% [{np.percentile(mbar, 2.5):+.4f}, '
              f'{np.percentile(mbar, 97.5):+.4f}],  worst D 95% [{np.percentile(mmax, 2.5):+.4f}, '
              f'{np.percentile(mmax, 97.5):+.4f}]   measured {MEAS["mu"]:+.3f}')
    # attribution for g: capacitor only vs everything
    S = simulate(5e-3, 5e-3, capTol=0.0)
    print('\n  g without capacitor tolerance, 95%:',
          np.round(np.percentile(100 * (S['g'] / 50 - 1), [2.5, 97.5]), 2))
    S = simulate(5e-3, 5e-3, resTol=0.0, scaleTol=0.0)
    print('  g with only capacitor tolerance, 95%:',
          np.round(np.percentile(100 * (S['g'] / 50 - 1), [2.5, 97.5]), 2))
    # rail sensitivity
    for vs in (4.9, 5.1):
        S = simulate(0, 0, VS=vs, resTol=0, capTol=0, scaleTol=0)
        print(f'  rail {vs} V (nominal parts):',
              {k: round(float(100 * (np.mean(S[k]) / DESIGN[k] - 1)), 2) for k in DESIGN},
              'mu mean', round(float(S['mu'].mean()), 4))
