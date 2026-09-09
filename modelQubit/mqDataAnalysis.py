# -*- coding: utf-8 -*-
"""
Plotting Model Qubit 
Benjamin Krawciw
8/26/25

The goal here is to take the data gathered from the model qubit and 
show via plots that it properly recreates the statistics of a qubit.
"""

#Library imports 
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#Setup variables
fname = 'modelQubit_09-05-2026_231938_noclip.log'

'''
TAKE IN THE DATA
'''

#Create lists to store the data
t, c, s, of, th, ph, x, y = ([] for i in range(8))

#Open the file
with open(fname, 'r') as f:
    #Create a reader object
    fRead = csv.reader(f)
    
    #Skip the headers
    next(fRead)
    next(fRead)
    
    #Record data
    for row in fRead:
        t.append(float(row[0]))
        c.append(float(row[1]))
        s.append(float(row[2]))
        #The extra float conversion before helps python ignore the space 
        of.append(bool(float((row[3]))))
        th.append(float(row[4]))
        ph.append(float(row[5]))
        x.append(float(row[6]))
        y.append(float(row[7]))

#Change lists to numpy arrays
t = np.array(t)
c = np.array(c)
s = np.array(s)
of = np.array(of)
th = np.array(th)
ph = np.array(ph)
x = np.array(x)
y = np.array(y)


'''
SORT THROUGH DATA, COMPUTE STATISTICS
'''

#Count how many distinct measurements were taken
mIndices = []
incycle = False
startdex = 0
#Iterate through on/off
for i, elem in enumerate(of):
    #Don't start looking until measurement activates
    if (c[i] == 0.0) and (s[i] == 0.0):
        continue

    #Measurement cycles begin and end when of switches on
    if (elem == True) and (of[i-1] == False):
        #For the first cycle,
        if incycle == False:
            startdex = i
            incycle = True
        #For subsequent cycles
        else:
            #End and record the previous cycle
            enddex = i-1
            mIndices.append( (startdex, enddex) )
            
            #Begin the next cycle
            startdex = i

#Close the final cycle at end-of-file (logs have variable length)
if incycle:
    mIndices.append( (startdex, len(of) - 1) )

count = len(mIndices)
print(f"Found {count} measurement cycles in {fname}")

#Normalize a sine-cosine pair
def normPair(cIn, sIn):
    #Ensure normalization
    csMag = np.sqrt((cIn)**2 + (sIn)**2)
    ci = cIn / csMag
    si = sIn / csMag
    return(ci, si)

#For each measurement, what was the relative angle between (x, y) and (c, s)?
def rangle(index):
    #Ensure normalization
    ci, si = normPair(c[index], s[index])
    xi, yi = normPair(x[index], y[index])
    
    #Find relative sine and cosine
    sr = yi * ci - xi * si
    cr = xi * ci + yi * si
    
    #Use arctan to find relative angle
    angle = np.arctan2(sr, cr)
    return(angle)

rangles = np.array([-rangle(dex[0]) for dex in mIndices])

#Did we measure up or down on this axis?
def updown(index, bound = 0.5):
    #Ensure normalization
    ci, si = normPair(c[index], s[index])
    xi, yi = normPair(x[index], y[index])
    
    #Find relative cosine
    cr = xi * ci + yi * si
    
    #Judge up or down within bound
    if cr > (1-bound):
        return 1
    elif cr < -(1-bound):
        return -1
    #Return 0 if inconclusive
    else:
        return 0

    
mmts = np.array([updown(i[1]) for i in mIndices])

#Bin the relative angles and compute measured probabilities
binCount = 10
binStep = (np.pi / 2.0) / (binCount - 1)
inBins = np.zeros(binCount) #Storing how many measurements in bin
bins = np.zeros(binCount) #Storing net measurement of bin
#Use every measurement found in this log (length varies by run)
nMmts = len(mIndices)
for i in range(nMmts):
    angle = rangles[i]
    #Which bin does the angle go in?
    binDex = int(np.round( angle / binStep ))
    #Skip angles that fall outside the [0, pi/2] analysis window
    if binDex < 0 or binDex >= binCount:
        continue
    #Count trial
    inBins[binDex] += 1
    #Count result
    bins[binDex] += (mmts[i] + 1.) / 2. 

#Compute probabilities of measuring up (empty bins -> nan)
probs = np.full(binCount, np.nan)
nonzero = inBins > 0
probs[nonzero] = bins[nonzero] / inBins[nonzero]
print(f"Using {nMmts} measurements")
print(inBins)
#print(probs)

#Compare to ideal probs
binThets = np.linspace(0.0, np.pi / 2., binCount)
ideal = np.cos(0.5 * binThets)**2.
#print(ideal)
#print(probs-ideal)

valid = ~np.isnan(probs)

'''
PLOT A SAMPLE OF DATA COLLECTED DURING TRIALS
'''
testPlot, testAx = plt.subplots(nrows = 2, ncols = 1)
nStart, nEnd = mIndices[3][0], mIndices[5][1]
testAx[0].plot(t[nStart:nEnd], c[nStart:nEnd] * of[nStart:nEnd], label = 'c')
testAx[0].plot(t[nStart:nEnd], s[nStart:nEnd] * of[nStart:nEnd], label = 's')
testAx[0].plot(t[nStart:nEnd], ph[nStart:nEnd], label = r'$\phi$')
testAx[1].plot(t[nStart:nEnd], x[nStart:nEnd], label = r'x', ls = '--')
testAx[1].plot(t[nStart:nEnd], y[nStart:nEnd], label = r'y', ls = '--')
testAx[1].plot(t[nStart:nEnd], th[nStart:nEnd], label = r'$\theta$', ls = '--')
testAx[0].legend()
testAx[1].legend()
testPlot.show()
testAx[1].set(xlabel = 't(s)',
              ylabel = 'Measurements')
testAx[0].set(ylabel = 'Control Signals')
testPlot.savefig('trialCycle.pdf')


'''
PLOT MEASURED VS IDEAL QUBIT
'''
#Creating the plot & axis objects
fig, ax = plt.subplots()

#Draw the plot of measured probabilities
ax.errorbar(binThets[valid], probs[valid], yerr = 0.1, marker = '.', linestyle = 'none',
            capsize = 5., ecolor = 'k', markersize = 20, color = 'red',
            label = 'empirical')

#Draw the ideal probabilities
thetaDom = np.linspace(0.0, np.pi / 2., 1000)
ideal = np.cos(0.5 * thetaDom)**2
ax.plot(thetaDom, ideal,
        label = 'ideal')

#Axes labels
ax.set_xlabel("Relative angle (rad)")
ax.set_ylabel("Up probability")

#Legend
ax.legend()

#Save figure
fig.savefig("empiricalProb.pdf")
fig.show()

#Compute and report reduced chi-squared of the experiment
#Only bins that received at least one count contribute (empty bins vary with log length)
dof = int(np.sum(valid)) #Note that no post-hoc fitting was done
var = 0.01 #Variance set by assumed error of 10%
theoreticals = np.cos(0.5 * binThets)**2
redChi = np.sum((probs[valid] - theoreticals[valid])**2) / (dof * var)
print("Reduced chi-squared: %0.2f" %redChi)


'''
ODE FIDELITY: does (x, y) collapse to sign(theta) * (c, s)?

Measurement force is the on/off column (of), not (c, s) leaving (0, 0).
Within each cycle, (c, s) stays fixed at the axis while of is True for the
first half (force on) and False for the second; (c, s) = (0, 0) only appears
at the start of the log before the first measurement.

Predictor: sign of raw logged theta at the first force-on sample (the analog
computer's theta voltage), not the mean over the collapse window.
'''

def wrapAngle(a):
    """Wrap angle(s) to (-pi, pi]."""
    return (a + np.pi) % (2.0 * np.pi) - np.pi

# Tolerance for "collapsed to the predicted pole"
angleTolDeg = 20.0
angleTol = np.deg2rad(angleTolDeg)

predAngles = []
finalAngles = []
initThetas = []
angleErrs = []

for start, end in mIndices:
    forceOn = of[start:end + 1]
    if not np.any(forceOn):
        continue

    onIdx = np.where(forceOn)[0]
    startIdx = start + onIdx[0]
    forceIdx = start + onIdx[-1]

    # Raw theta at the start of the collapse (analog theta voltage)
    th0 = th[startIdx]
    initThetas.append(th0)

    # Final state: last sample while measurement force is on
    finalAng = np.arctan2(y[forceIdx], x[forceIdx])
    axisAng = np.arctan2(s[forceIdx], c[forceIdx])

    # Positive initial theta -> +(c, s); negative -> -(c, s)
    if th0 >= 0.0:
        predAng = axisAng
    else:
        predAng = wrapAngle(axisAng + np.pi)

    predAngles.append(predAng)
    finalAngles.append(finalAng)
    angleErrs.append(wrapAngle(finalAng - predAng))

predAngles = np.array(predAngles)
finalAngles = np.array(finalAngles)
initThetas = np.array(initThetas)
angleErrs = np.array(angleErrs)

success = np.abs(angleErrs) <= angleTol
nFid = len(success)
nOk = int(np.sum(success))
successPct = 100.0 * nOk / nFid if nFid else float('nan')

print("\n--- ODE collapse fidelity ---")
ofFrac = float(np.mean([of[a:b + 1].mean() for a, b in mIndices]))
print(f"Force-on window: of == True ({100.0 * ofFrac:.0f}% of each cycle on average)")
print(f"Predictor: sign of raw theta at force-on start")
print(f"|theta_0|: success={np.mean(np.abs(initThetas[success])):.3f}, "
      f"fail={np.mean(np.abs(initThetas[~success])):.3f}")
print(f"Collapse success (within +/- {angleTolDeg:.0f} deg of sign(theta_0)*(c,s)): "
      f"{nOk}/{nFid} = {successPct:.1f}%")
print(f"  among theta_0 >= 0: {int(np.sum(success[initThetas >= 0]))}/"
      f"{int(np.sum(initThetas >= 0))}")
print(f"  among theta_0 <  0: {int(np.sum(success[initThetas < 0]))}/"
      f"{int(np.sum(initThetas < 0))}")

# Scatter: predicted pole angle vs final (x, y) angle
fidFig, fidAx = plt.subplots()
fidAx.scatter(predAngles[success], finalAngles[success],
              c = 'tab:green', marker = '.', label = 'success', zorder = 2)
fidAx.scatter(predAngles[~success], finalAngles[~success],
              c = 'tab:red', marker = '.', label = 'fail', zorder = 2)
# Ideal: final == predicted (mod 2pi wrapping makes a single diagonal imperfect)
diag = np.linspace(-np.pi, np.pi, 200)
fidAx.plot(diag, diag, 'k--', lw = 1, label = r'final = predicted', zorder = 1)
fidAx.set_xlabel(r'Predicted angle of $\mathrm{sign}(\theta_0)\,(c,s)$ (rad)')
fidAx.set_ylabel(r'Final angle of $(x,y)$ (rad)')
fidAx.set_xlim(-np.pi, np.pi)
fidAx.set_ylim(-np.pi, np.pi)
fidAx.set_aspect('equal')
fidAx.legend()
fidFig.tight_layout()
fidFig.savefig('odeFidelityScatter.pdf')
fidFig.show()


'''
ODE RESIDUALS VIA FINITE DIFFERENCES

Model (analog computer), in natural units where the circle is r = 1:
  dx/dt = v_x,  dy/dt = v_y
  dv_x/dt = k x (1 - r^2) - d v_x + theta * c
  dv_y/dt = k y (1 - r^2) - d v_y + theta * s

Logged x, y, c, s are ADC fractions of 5 V with observed |r| ≈ 0.87, not 1, so
(1 - r_ADC^2) does not vanish on the physical circle. Rescale
  (x, y, c, s)_nat = (x, y, c, s)_ADC / R_obs
with R_obs = median |r|_ADC, then fit in natural units:
  a ≈ k q (1 - r_nat^2) - d v + g theta n_eff
Theta keeps its own voltage scale; g absorbs that mismatch.

v is not logged: estimate v and a = dv/dt from finite differences in t.
Firmware: when of is False, measureOff() drives c/s pins to 0 V, but the log
still stores the software axis. Use n_eff = of * (c, s)_nat.
'''

# Observed circle radius in logged (ADC fraction) units
rAdc = np.sqrt(x**2 + y**2)
Robs = float(np.median(rAdc))
xN, yN = x / Robs, y / Robs
cN, sN = c / Robs, s / Robs

# Effective measurement axis on the wires (natural units)
cEff = cN * of.astype(float)
sEff = sN * of.astype(float)

# Finite-difference kinematics in natural units
vx = np.gradient(xN, t)
vy = np.gradient(yN, t)
ax = np.gradient(vx, t)
ay = np.gradient(vy, t)

r2 = xN**2 + yN**2
# Drop endpoints: gradient one-sided there is less trustworthy
interior = np.zeros(len(t), dtype = bool)
interior[1:-1] = True

def odeDesign(q, v, nEff, thIn, r2In, mask):
    """Columns for [k, d, g] in a ≈ k q(1-r^2) - d v + g theta n."""
    qq = q[mask]
    return np.column_stack((
        qq * (1.0 - r2In[mask]),
        -v[mask],
        thIn[mask] * nEff[mask],
    ))

# Stack x and y equations on interior samples
mask = interior
Ax = odeDesign(xN, vx, cEff, th, r2, mask)
Ay = odeDesign(yN, vy, sEff, th, r2, mask)
Amat = np.vstack((Ax, Ay))
bvec = np.concatenate((ax[mask], ay[mask]))

coeff, _, _, _ = np.linalg.lstsq(Amat, bvec, rcond = None)
kFit, dFit, gFit = coeff
aPred = Amat @ coeff
resid = bvec - aPred

def rSquared(yTrue, yHat):
    ssRes = np.sum((yTrue - yHat)**2)
    ssTot = np.sum((yTrue - np.mean(yTrue))**2)
    return 1.0 - ssRes / ssTot if ssTot > 0 else float('nan')

r2All = rSquared(bvec, aPred)

# Same fit restricted to force-on samples only
maskOn = interior & of
AxOn = odeDesign(xN, vx, cEff, th, r2, maskOn)
AyOn = odeDesign(yN, vy, sEff, th, r2, maskOn)
Aon = np.vstack((AxOn, AyOn))
bon = np.concatenate((ax[maskOn], ay[maskOn]))
coeffOn, _, _, _ = np.linalg.lstsq(Aon, bon, rcond = None)
kOn, dOn, gOn = coeffOn
aPredOn = Aon @ coeffOn
r2On = rSquared(bon, aPredOn)

# Free radial target in natural units: a ≈ b1 q - b2 q r^2 - d v + g theta n
# For k q (R^2 - r^2): b1 = k R^2, b2 = k, so R^2 = b1/b2 (expect ~1 after rescaling)
def odeDesignFreeR(q, v, nEff, thIn, r2In, mask):
    qq = q[mask]
    return np.column_stack((
        qq,
        -qq * r2In[mask],
        -v[mask],
        thIn[mask] * nEff[mask],
    ))

Afree = np.vstack((
    odeDesignFreeR(xN, vx, cEff, th, r2, mask),
    odeDesignFreeR(yN, vy, sEff, th, r2, mask),
))
coeffFree, _, _, _ = np.linalg.lstsq(Afree, bvec, rcond = None)
b1, b2, dFree, gFree = coeffFree
r2Target = b1 / b2 if abs(b2) > 1e-12 else float('nan')
aPredFree = Afree @ coeffFree
r2Free = rSquared(bvec, aPredFree)

print("\n--- ODE finite-difference residuals ---")
print(f"Rescaled to natural units by R_obs = median |r|_ADC = {Robs:.3f}")
print(f"Samples used (interior): {int(np.sum(mask))} "
      f"({int(np.sum(maskOn))} force-on); median dt = {np.median(np.diff(t)):.3f} s")
print(f"Fit a ≈ k q(1-r_nat^2) - d v + g theta n_eff  [all interior]:")
print(f"  k = {kFit:.3g},  d = {dFit:.3g},  g = {gFit:.3g},  R^2 = {r2All:.3f}")
print(f"  (k>0 => radial restoring toward unit circle; d>0 => drag)")
print(f"Fit same form [force-on only]:")
print(f"  k = {kOn:.3g},  d = {dOn:.3g},  g = {gOn:.3g},  R^2 = {r2On:.3f}")
print(f"Fit a ≈ b1 q - b2 q r^2 - d v + g theta n_eff  [free radial target]:")
print(f"  b1 = {b1:.3g},  b2 = {b2:.3g},  d = {dFree:.3g},  g = {gFree:.3g}")
print(f"  implied R_target^2 = b1/b2 = {r2Target:.3f}  (|R|_nat = {np.sqrt(abs(r2Target)):.3f}),  "
      f"R^2 = {r2Free:.3f}")
print(f"Residual std (all-interior unit-circle fit): {np.std(resid):.3g}")

# Predicted vs measured acceleration scatter
odeFig, odeAx = plt.subplots()
# Subsample for readability
rng = np.random.default_rng(0)
nShow = min(4000, len(bvec))
show = rng.choice(len(bvec), size = nShow, replace = False)
odeAx.scatter(bvec[show], aPred[show], s = 6, alpha = 0.25, label = 'samples')
lim = np.percentile(np.abs(bvec), 99)
odeAx.plot([-lim, lim], [-lim, lim], 'k--', lw = 1, label = 'predicted = measured')
odeAx.set_xlabel(r'Measured $a$ from finite differences')
odeAx.set_ylabel(r'Predicted $a = k q(1-r^2) - d v + g\theta n$')
odeAx.set_title(r'ODE acceleration check, raw ($R^2 = %0.3f$)' % r2All)
odeAx.set_aspect('equal', adjustable = 'box')
odeAx.set_xlim(-lim, lim)
odeAx.set_ylim(-lim, lim)
odeAx.legend()
odeFig.tight_layout()
odeFig.savefig('odeResidualScatter.pdf')
odeFig.show()

# Time-series residual snapshot over a few cycles
resFig, resAx = plt.subplots(nrows = 2, ncols = 1, sharex = True)
# Reconstruct per-sample prediction on interior with the all-interior fit
aPredX = np.full(len(t), np.nan)
aPredY = np.full(len(t), np.nan)
aPredX[mask] = Ax @ coeff
aPredY[mask] = Ay @ coeff
n0, n1 = mIndices[3][0], mIndices[5][1]
resAx[0].plot(t[n0:n1], ax[n0:n1], label = r'$a_x$ meas')
resAx[0].plot(t[n0:n1], aPredX[n0:n1], label = r'$a_x$ ODE', ls = '--')
resAx[1].plot(t[n0:n1], ay[n0:n1], label = r'$a_y$ meas')
resAx[1].plot(t[n0:n1], aPredY[n0:n1], label = r'$a_y$ ODE', ls = '--')
resAx[0].legend(); resAx[1].legend()
resAx[1].set_xlabel('t (s)')
resAx[0].set_ylabel(r'$a_x$')
resAx[1].set_ylabel(r'$a_y$')
resFig.suptitle('Raw finite differences')
resFig.tight_layout()
resFig.savefig('odeResidualTimeseries.pdf')
resFig.show()


'''
ODE RESIDUALS WITH SMOOTHED TRAJECTORIES

Smooth natural-unit x, y and raw theta with a centered moving average, then
re-estimate v and a. RHS uses the same smoothed state (and of*axis).
'''

def movingAverage(arr, window):
    """Centered moving average; window must be odd. Edges use partial windows."""
    if window % 2 == 0:
        raise ValueError('window must be odd')
    w = window
    # np.convolve with 'same' uses zero-padding; correct with a count kernel
    num = np.convolve(arr, np.ones(w), mode = 'same')
    den = np.convolve(np.ones_like(arr, dtype = float), np.ones(w), mode = 'same')
    return num / den

def fitOdeAcceleration(qx, qy, thIn, vxIn, vyIn, axIn, ayIn, r2In, valid):
    """Least-squares fit of a ≈ k q(1-r^2) - d v + g theta n_eff on valid samples."""
    maskFit = valid
    maskFitOn = valid & of
    Aall = np.vstack((
        odeDesign(qx, vxIn, cEff, thIn, r2In, maskFit),
        odeDesign(qy, vyIn, sEff, thIn, r2In, maskFit),
    ))
    ball = np.concatenate((axIn[maskFit], ayIn[maskFit]))
    coeffAll, _, _, _ = np.linalg.lstsq(Aall, ball, rcond = None)
    predAll = Aall @ coeffAll

    Aon = np.vstack((
        odeDesign(qx, vxIn, cEff, thIn, r2In, maskFitOn),
        odeDesign(qy, vyIn, sEff, thIn, r2In, maskFitOn),
    ))
    bon = np.concatenate((axIn[maskFitOn], ayIn[maskFitOn]))
    coeffOnLoc, _, _, _ = np.linalg.lstsq(Aon, bon, rcond = None)
    predOn = Aon @ coeffOnLoc

    Afr = np.vstack((
        odeDesignFreeR(qx, vxIn, cEff, thIn, r2In, maskFit),
        odeDesignFreeR(qy, vyIn, sEff, thIn, r2In, maskFit),
    ))
    coeffFr, _, _, _ = np.linalg.lstsq(Afr, ball, rcond = None)
    predFr = Afr @ coeffFr
    b1Loc, b2Loc, dLoc, gLoc = coeffFr
    r2T = b1Loc / b2Loc if abs(b2Loc) > 1e-12 else float('nan')

    return {
        'coeff': coeffAll,
        'r2': rSquared(ball, predAll),
        'resid_std': float(np.std(ball - predAll)),
        'coeff_on': coeffOnLoc,
        'r2_on': rSquared(bon, predOn),
        'coeff_free': coeffFr,
        'r2_free': rSquared(ball, predFr),
        'r2_target': r2T,
        'pred': predAll,
        'meas': ball,
        'A_x': odeDesign(qx, vxIn, cEff, thIn, r2In, maskFit),
        'A_y': odeDesign(qy, vyIn, sEff, thIn, r2In, maskFit),
        'mask': maskFit,
    }

smoothWindows = [5, 7, 11, 15]
print("\n--- ODE residuals with smoothed trajectories (natural units) ---")
smoothResults = {}
for win in smoothWindows:
    xs = movingAverage(xN, win)
    ys = movingAverage(yN, win)
    ths = movingAverage(th, win)
    r2s = xs**2 + ys**2
    vxs = np.gradient(xs, t)
    vys = np.gradient(ys, t)
    axs = np.gradient(vxs, t)
    ays = np.gradient(vys, t)
    # Exclude edges heavily influenced by partial smoothing windows
    half = win // 2
    valid = np.zeros(len(t), dtype = bool)
    valid[half + 1:-(half + 1)] = True
    out = fitOdeAcceleration(xs, ys, ths, vxs, vys, axs, ays, r2s, valid)
    smoothResults[win] = dict(out, xs = xs, ys = ys, ths = ths,
                              vxs = vxs, vys = vys, axs = axs, ays = ays)
    kS, dS, gS = out['coeff']
    print(f"window = {win} samples (~{win * np.median(np.diff(t)):.2f} s):")
    print(f"  unit-circle fit: k={kS:.3g}, d={dS:.3g}, g={gS:.3g}, "
          f"R^2={out['r2']:.3f}, resid_std={out['resid_std']:.3g}")
    print(f"  force-on only:   k={out['coeff_on'][0]:.3g}, d={out['coeff_on'][1]:.3g}, "
          f"g={out['coeff_on'][2]:.3g}, R^2={out['r2_on']:.3f}")
    print(f"  free radial:     |R|_nat={np.sqrt(abs(out['r2_target'])):.3f}, "
          f"R^2={out['r2_free']:.3f}")

# Use the window with best unit-circle R^2 for plots
bestWin = max(smoothResults, key = lambda w: smoothResults[w]['r2'])
best = smoothResults[bestWin]
print(f"Best smoothing window by R^2: {bestWin}")

odeSmFig, odeSmAx = plt.subplots()
measSm, predSm = best['meas'], best['pred']
nShowSm = min(4000, len(measSm))
showSm = rng.choice(len(measSm), size = nShowSm, replace = False)
odeSmAx.scatter(measSm[showSm], predSm[showSm], s = 6, alpha = 0.25, label = 'samples')
limSm = np.percentile(np.abs(measSm), 99)
odeSmAx.plot([-limSm, limSm], [-limSm, limSm], 'k--', lw = 1, label = 'predicted = measured')
odeSmAx.set_xlabel(r'Measured $a$ (smoothed finite differences)')
odeSmAx.set_ylabel(r'Predicted $a = k q(1-r^2) - d v + g\theta n$')
odeSmAx.set_title(r'Smoothed ODE check, window=%d ($R^2 = %0.3f$)' % (bestWin, best['r2']))
odeSmAx.set_aspect('equal', adjustable = 'box')
odeSmAx.set_xlim(-limSm, limSm)
odeSmAx.set_ylim(-limSm, limSm)
odeSmAx.legend()
odeSmFig.tight_layout()
odeSmFig.savefig('odeResidualScatterSmooth.pdf')
odeSmFig.show()

resSmFig, resSmAx = plt.subplots(nrows = 2, ncols = 1, sharex = True)
aPredXsm = np.full(len(t), np.nan)
aPredYsm = np.full(len(t), np.nan)
aPredXsm[best['mask']] = best['A_x'] @ best['coeff']
aPredYsm[best['mask']] = best['A_y'] @ best['coeff']
resSmAx[0].plot(t[n0:n1], best['axs'][n0:n1], label = r'$a_x$ smooth meas')
resSmAx[0].plot(t[n0:n1], aPredXsm[n0:n1], label = r'$a_x$ ODE', ls = '--')
resSmAx[1].plot(t[n0:n1], best['ays'][n0:n1], label = r'$a_y$ smooth meas')
resSmAx[1].plot(t[n0:n1], aPredYsm[n0:n1], label = r'$a_y$ ODE', ls = '--')
resSmAx[0].legend(); resSmAx[1].legend()
resSmAx[1].set_xlabel('t (s)')
resSmAx[0].set_ylabel(r'$a_x$')
resSmAx[1].set_ylabel(r'$a_y$')
resSmFig.suptitle(f'Smoothed finite differences (window={bestWin})')
resSmFig.tight_layout()
resSmFig.savefig('odeResidualTimeseriesSmooth.pdf')
resSmFig.show()


'''
ODE CHECKS THAT AVOID POINTWISE ACCELERATION

Acceleration residuals amplify timestep noise. Keep the a-based fits above for
scale / (k, d, g) exploration, but score ODE obedience with:

1) Integrated velocity increments (one derivative only):
     v(t+Tau) - v(t)  ≈  ∫_t^{t+Tau} [k q(1-r^2) - d v + g theta n] dt
   Fit (k, d, g) by least squares on these increments (trapezoid integrals).

2) Short-horizon trajectory prediction: from each force-on start, integrate
     dq/dt = v,  dv/dt = k q(1-r^2) - d v + g theta n_eff
   with measured theta(t) and n_eff(t) as inputs; compare to measured (x, y).
'''

def trapzIntegrate(ySeg, tSeg):
    """Trapezoidal integral of ySeg over tSeg (1-D)."""
    ySeg = np.asarray(ySeg, dtype = float)
    tSeg = np.asarray(tSeg, dtype = float)
    return float(np.sum(0.5 * (ySeg[1:] + ySeg[:-1]) * np.diff(tSeg)))

def fitIntegratedDeltaV(qx, qy, thIn, vxIn, vyIn, r2In, valid, lag, label):
    """
    Fit k, d, g from Delta v ≈ Int RHS over `lag` samples.
    Returns coeff, R^2, and arrays of measured/predicted Delta v.
    """
    starts = np.where(valid[:-lag] & valid[lag:])[0]
    # Prefer intervals that stay force-on or stay force-off (avoid of edges)
    rows = []
    dvMeas = []
    for i in starts:
        j = i + lag
        # require of constant on the open interval of samples used
        ofSeg = of[i:j + 1]
        if not (np.all(ofSeg) or np.all(~ofSeg)):
            continue
        tSeg = t[i:j + 1]
        for q, vcomp, nEff in (
            (qx, vxIn, cEff),
            (qy, vyIn, sEff),
        ):
            intRad = trapzIntegrate(q[i:j + 1] * (1.0 - r2In[i:j + 1]), tSeg)
            intDrag = trapzIntegrate(-vcomp[i:j + 1], tSeg)
            intForce = trapzIntegrate(thIn[i:j + 1] * nEff[i:j + 1], tSeg)
            rows.append((intRad, intDrag, intForce))
            dvMeas.append(vcomp[j] - vcomp[i])

    if len(rows) < 10:
        print(f"{label}: too few intervals (lag={lag})")
        return None

    Aint = np.array(rows, dtype = float)
    bint = np.array(dvMeas, dtype = float)
    coeffInt, _, _, _ = np.linalg.lstsq(Aint, bint, rcond = None)
    predInt = Aint @ coeffInt
    r2Int = rSquared(bint, predInt)
    kI, dI, gI = coeffInt
    print(f"{label} (lag={lag} samples, ~{lag * np.median(np.diff(t)):.2f} s, "
          f"n={len(bint)} increments):")
    print(f"  k={kI:.3g}, d={dI:.3g}, g={gI:.3g}, R^2={r2Int:.3f}")
    return {
        'coeff': coeffInt,
        'r2': r2Int,
        'meas': bint,
        'pred': predInt,
        'lag': lag,
    }

print("\n--- ODE integrated Delta-v checks (natural units) ---")
# Use lightly smoothed kinematics for v (still one derivative, much stabler than a)
softWin = 5
xSoft = movingAverage(xN, softWin)
ySoft = movingAverage(yN, softWin)
thSoft = movingAverage(th, softWin)
r2Soft = xSoft**2 + ySoft**2
vxSoft = np.gradient(xSoft, t)
vySoft = np.gradient(ySoft, t)
halfSoft = softWin // 2
validSoft = np.zeros(len(t), dtype = bool)
validSoft[halfSoft + 1:-(halfSoft + 1)] = True

intFits = {}
for lag in (3, 5, 8, 10):
    fit = fitIntegratedDeltaV(
        xSoft, ySoft, thSoft, vxSoft, vySoft, r2Soft, validSoft, lag,
        label = f"Integrated Delta-v",
    )
    if fit is not None:
        intFits[lag] = fit

bestLag = max(intFits, key = lambda L: intFits[L]['r2']) if intFits else None
if bestLag is not None:
    bestInt = intFits[bestLag]
    print(f"Best integrated fit by R^2: lag={bestLag}")
    kTraj, dTraj, gTraj = bestInt['coeff']

    dvFig, dvAx = plt.subplots()
    nDv = min(4000, len(bestInt['meas']))
    showDv = rng.choice(len(bestInt['meas']), size = nDv, replace = False)
    dvAx.scatter(bestInt['meas'][showDv], bestInt['pred'][showDv],
                 s = 6, alpha = 0.25, label = 'increments')
    limDv = np.percentile(np.abs(bestInt['meas']), 99)
    dvAx.plot([-limDv, limDv], [-limDv, limDv], 'k--', lw = 1,
              label = 'predicted = measured')
    dvAx.set_xlabel(r'Measured $\Delta v$')
    dvAx.set_ylabel(r'Predicted $\int (k q(1-r^2)-d v+g\theta n)\,dt$')
    dvAx.set_title(r'Integrated $\Delta v$ check, lag=%d ($R^2=%0.3f$)' %
                   (bestLag, bestInt['r2']))
    dvAx.set_aspect('equal', adjustable = 'box')
    dvAx.set_xlim(-limDv, limDv)
    dvAx.set_ylim(-limDv, limDv)
    dvAx.legend()
    dvFig.tight_layout()
    dvFig.savefig('odeDeltaVScatter.pdf')
    dvFig.show()
else:
    kTraj, dTraj, gTraj = kOn, dOn, gOn
    print("Falling back to force-on acceleration-fit coefficients for trajectories")


def integratePlanar(i0, i1, kP, dP, gP):
    """
    Integrate 2D ODE from sample i0 to i1 using measured theta and n_eff.
    Initial (q, v) from smoothed natural-unit kinematics at i0.
    Returns predicted x,y arrays aligned with t[i0:i1+1].
    """
    nSteps = i1 - i0
    qx = np.empty(nSteps + 1)
    qy = np.empty(nSteps + 1)
    vxP = np.empty(nSteps + 1)
    vyP = np.empty(nSteps + 1)
    qx[0], qy[0] = xSoft[i0], ySoft[i0]
    vxP[0], vyP[0] = vxSoft[i0], vySoft[i0]

    for s in range(nSteps):
        i = i0 + s
        dtStep = t[i + 1] - t[i]
        th0, th1 = thSoft[i], thSoft[i + 1]
        cx0, cx1 = cEff[i], cEff[i + 1]
        sy0, sy1 = sEff[i], sEff[i + 1]

        def accel(qx_, qy_, vx_, vy_, th_, cx_, sy_):
            r2_ = qx_ * qx_ + qy_ * qy_
            ax_ = kP * qx_ * (1.0 - r2_) - dP * vx_ + gP * th_ * cx_
            ay_ = kP * qy_ * (1.0 - r2_) - dP * vy_ + gP * th_ * sy_
            return ax_, ay_

        # RK2 with inputs averaged over the step
        thM = 0.5 * (th0 + th1)
        cxM = 0.5 * (cx0 + cx1)
        syM = 0.5 * (sy0 + sy1)
        ax1, ay1 = accel(qx[s], qy[s], vxP[s], vyP[s], th0, cx0, sy0)
        qxMid = qx[s] + 0.5 * dtStep * vxP[s]
        qyMid = qy[s] + 0.5 * dtStep * vyP[s]
        vxMid = vxP[s] + 0.5 * dtStep * ax1
        vyMid = vyP[s] + 0.5 * dtStep * ay1
        ax2, ay2 = accel(qxMid, qyMid, vxMid, vyMid, thM, cxM, syM)
        qx[s + 1] = qx[s] + dtStep * vxMid
        qy[s + 1] = qy[s] + dtStep * vyMid
        vxP[s + 1] = vxP[s] + dtStep * ax2
        vyP[s + 1] = vyP[s] + dtStep * ay2

    return qx, qy


print("\n--- ODE trajectory prediction (force-on windows) ---")
trajRmse = []
angErrDeg = []
finalPosErr = []
trajExamples = []  # a few (t, xy_meas, xy_pred) for plotting

for ci, (start, end) in enumerate(mIndices):
    forceOn = of[start:end + 1]
    if not np.any(forceOn):
        continue
    onIdx = np.where(forceOn)[0]
    i0 = start + onIdx[0]
    i1 = start + onIdx[-1]
    if i1 - i0 < 3:
        continue
    if not (validSoft[i0] and validSoft[i1]):
        continue

    qxP, qyP = integratePlanar(i0, i1, kTraj, dTraj, gTraj)
    qxM, qyM = xSoft[i0:i1 + 1], ySoft[i0:i1 + 1]
    err = np.sqrt((qxP - qxM)**2 + (qyP - qyM)**2)
    trajRmse.append(float(np.mean(err)))
    finalPosErr.append(float(err[-1]))

    angM = np.arctan2(qyM[-1], qxM[-1])
    angP = np.arctan2(qyP[-1], qxP[-1])
    angErrDeg.append(float(np.degrees(wrapAngle(angP - angM))))

    if ci in (3, 8, 89, 94) or (len(trajExamples) < 1):
        trajExamples.append((ci, t[i0:i1 + 1], qxM, qyM, qxP, qyP))

trajRmse = np.array(trajRmse)
angErrDeg = np.array(angErrDeg)
finalPosErr = np.array(finalPosErr)
print(f"Using k,d,g from integrated fit (lag={bestLag}): "
      f"k={kTraj:.3g}, d={dTraj:.3g}, g={gTraj:.3g}")
print(f"Force-on trajectories scored: {len(trajRmse)}")
print(f"Path RMSE (natural units): mean={np.mean(trajRmse):.3f}, "
      f"median={np.median(trajRmse):.3f}")
print(f"Final |Δr| (natural units): mean={np.mean(finalPosErr):.3f}, "
      f"median={np.median(finalPosErr):.3f}")
print(f"Final angle error |Δψ|: mean={np.mean(np.abs(angErrDeg)):.1f} deg, "
      f"median={np.median(np.abs(angErrDeg)):.1f} deg")
print(f"Final angle within +/- 20 deg: "
      f"{np.mean(np.abs(angErrDeg) <= 20.0) * 100:.1f}%")

# Histogram of final angle errors
angFig, angAx = plt.subplots()
angAx.hist(angErrDeg, bins = 36, range = (-180, 180), color = 'tab:blue', alpha = 0.8)
angAx.axvline(0.0, color = 'k', ls = '--', lw = 1)
angAx.set_xlabel(r'Final angle error $\psi_{\mathrm{pred}}-\psi_{\mathrm{meas}}$ (deg)')
angAx.set_ylabel('Counts')
angAx.set_title('ODE trajectory prediction (force-on)')
angFig.tight_layout()
angFig.savefig('odeTrajAngleError.pdf')
angFig.show()

# Example overlays in the plane
exFig, exAxes = plt.subplots(nrows = 2, ncols = 2, figsize = (7, 7),
                             subplot_kw = {'aspect': 'equal'})
exAxes = exAxes.ravel()
unitCirc = np.linspace(0, 2 * np.pi, 200)
for axE, ex in zip(exAxes, trajExamples[:4]):
    ci, tt, qxM, qyM, qxP, qyP = ex
    axE.plot(np.cos(unitCirc), np.sin(unitCirc), 'k-', lw = 0.5, alpha = 0.4)
    axE.plot(qxM, qyM, label = 'meas', lw = 2)
    axE.plot(qxP, qyP, label = 'ODE', ls = '--', lw = 2)
    axE.plot(qxM[0], qyM[0], 'o', color = 'k', markersize = 4)
    axE.set_xlim(-1.5, 1.5); axE.set_ylim(-1.5, 1.5)
    axE.set_title(f'trial {ci}')
    axE.legend(fontsize = 8)
exFig.suptitle('Force-on trajectories: measured vs ODE integrate')
exFig.tight_layout()
exFig.savefig('odeTrajExamples.pdf')
exFig.show()


'''
ANIMATIONS
'''

#Plotting csxy state
fig, ax = plt.subplots(figsize = [3., 3.])
# Move the left and bottom spines to x = 0 and y = 0, respectively.
ax.spines[["left", "bottom"]].set_position(("data", 0))
# Hide the top and right spines.
ax.spines[["top", "right"]].set_visible(False)

# Draw arrows (as black triangles: ">k"/"^k") at the end of the axes.  In each
# case, one of the coordinates (0) is a data coordinate (i.e., y = 0 or x = 0,
# respectively) and the other one (1) is an axes coordinate (i.e., at the very
# right/top of the axes).  Also, disable clipping (clip_on=False) as the marker
# actually spills out of the Axes.
ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)
ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)
#Hide axis numbers, draw grid
ax.tick_params(grid_alpha = 0.5, 
               labelcolor = 'none',
               labelsize = 0.)

#Draw the unit circle
unitC = plt.Circle((0., 0.), 1., fill = False, ec = 'k')
ax.add_patch(unitC)

#Safe initial frame (old logs assumed index 1000 always existed)
i0 = min(1000, len(x) - 1)

#Plot an xy point
xy, = ax.plot([x[i0]], [y[i0]], marker = 'o')

#Draw the measurement axis
axup, = ax.plot([0., c[i0]], [0., s[i0]], color = 'orange')
axdown, = ax.plot([0., -c[i0]], [0., -s[i0]], color = 'orange', ls = '--')
axtip, = ax.plot(c[i0], s[i0], color = 'orange', marker = 'o', mfc = 'white')

#Indicator for measuring
ind = plt.Circle((0.95, 0.95), 0.05, color = 'red')
ax.add_patch(ind)

#Function for animating a frame
def step(i):
    #Move the xy point to the proper location
    xy.set_data(
        [x[i]], 
        [y[i]]
        )
    
    #Move the measurement axis to the proper location
    axup.set_data(
        [0., c[i]], 
        [0., s[i]]
        )
    axdown.set_data(
        [0., -c[i]],
        [0., -s[i]]
        )
    axtip.set_data(
        [c[i]],
        [s[i]]
        )

    #Turn the indicator on/off
    ind.set(visible = of[i])
    
    return (xy, axup, axdown, axtip, ind)

#Keep animation frames inside the recorded sample range
def trial_frames(start, end, pad = 40):
    lo = max(0, start - pad)
    hi = min(len(t) - 1, end + pad)
    return range(lo, hi + 1)

#Save / skip GIFs in animations/
animDir = 'animations'
os.makedirs(animDir, exist_ok = True)

#Animation of full trial
for i in range(len(mIndices)):
    ffname = os.path.join(animDir, 'test_' + str(i) + '.gif')
    if os.path.exists(ffname):
        print('Skipping existing video for trial ' + str(i))
        continue

    print('Making video for trial ' + str(i))
    
    start = mIndices[i][0]
    end = mIndices[i][1]
    video = ani.FuncAnimation(fig, step,
                              repeat = True, 
                              frames = trial_frames(start, end),
                              interval = 50,
                              blit = True)
    #fig.show()
    
    #To save the animation using Pillow as a gif
    writer = ani.PillowWriter(fps=20)
    video.save(ffname, writer=writer)

#Optional: animate a specific inconclusive trial if it exists in this log
#It's clear that it was on its way to the expected state, but the measurement
# period ended too early.
'''
weirdDex = 740
weirdPath = os.path.join(animDir, 'weirdTest5.gif')
if weirdDex < len(mIndices) and not os.path.exists(weirdPath):
    start = mIndices[weirdDex][0]
    end = mIndices[weirdDex][1]
    video = ani.FuncAnimation(fig, step,
                              repeat = True, 
                              frames = trial_frames(start, end),
                              interval = 50,
                              blit = True)
    writer = ani.PillowWriter(fps=20)
    video.save(weirdPath, writer=writer)
'''
