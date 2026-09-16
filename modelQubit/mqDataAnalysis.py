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
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#Setup variables
fname = 'modelQubit_09-15-2026_165442.log'

'''
TAKE IN THE DATA
'''

#Create lists to store the data
t, c, s, of, th, ph, x, y = ([] for i in range(8))

#Open the file
with open(fname, 'r') as f:
    #Create a reader object
    fRead = csv.reader(f)
    
    #Read the headers (first line identifies the firmware version)
    headerLine = ','.join(next(fRead))
    next(fRead)
    fwMatch = re.search(r'fw v(\d+)', headerLine)
    fwVer = int(fwMatch.group(1)) if fwMatch else 1
    fwV2 = fwVer >= 2
    
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
FIRMWARE-VERSION CORRECTIONS

Old firmware (before fw v2):
  - converted theta with 2*raw - 1 (the x/y level shifter), but the theta
    interface is V_pin = 2V + 0.8 V_theta. Correct: theta = V_theta/5V
    = 1.25*raw - 0.5 = 0.625*theta_old + 0.125.
  - printed 2 decimals, so logged values move in 0.01 steps.
  - read immediately after switching c/s, before the 10 ms PWM filter
    settled, so the first force-on theta is still the force-off theta.
fw v2 fixes all three in firmware (4 decimals, settle delay).
'''
if fwV2:
    printStep = 1e-4
    thetaLag = 0  # first force-on sample is already settled
else:
    th = 0.625 * th + 0.125
    printStep = 0.01
    thetaLag = 2  # skip the unsettled samples when reading theta_0
# Logged position resolution (fraction of 5V): coarser of ADC and print steps
logStep = max(2.0 / 1023.0, printStep)
print(f"Firmware: {f'v{fwVer}' if fwV2 else 'pre-v2 (theta converted, 2-decimal log)'}; "
      f"logged position step = {logStep:.4g}")

# measureOff() writes PWM 128, which is +1/255 of full scale above 0 V
C_OFF = 2.0 * 128.0 / 255.0 - 1.0


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
THETA CALIBRATION AND BORN-RULE PREDICTION

Fit the theta circuit's actual transfer on settled force-on samples:
  theta = a0 + b u + g phi,   u = x c + y s   (all in fractions of 5 V)
The ideal circuit has theta ∝ 1 + cos(D) - phi_nat, so with a unit axis and
the state at radius R (the 14/19 artifact), theta changes sign at
  phi = (a0 + b R cos D) / |g|.
A uniform phi sweep reproduces the Born rule exactly when it spans
  [ (a0 - b R)/|g| , (a0 + b R)/|g| ],
which is what the firmware's V_PHI_LO / V_PHI_HI should be set to.
'''
print("\n--- Theta calibration ---")
calIdx = []
for a_, b_ in mIndices:
    onIdx = a_ + np.where(of[a_:b_ + 1])[0]
    calIdx += list(onIdx if fwV2 else onIdx[1:])  # skip unsettled sample on old logs
calIdx = np.array(calIdx, dtype = int)
uCal = x[calIdx] * c[calIdx] + y[calIdx] * s[calIdx]
Acal = np.column_stack((np.ones(len(calIdx)), uCal, ph[calIdx]))
calCo, _, _, _ = np.linalg.lstsq(Acal, th[calIdx], rcond = None)
a0Cal, bCal, gCal = calCo
calRes = th[calIdx] - Acal @ calCo
calR2 = 1.0 - np.var(calRes) / np.var(th[calIdx])
rMed = float(np.median(np.hypot(x, y)))
axisMag = float(np.median(np.hypot(c[calIdx], s[calIdx])))
print(f"theta = {a0Cal:.4f} + {bCal:.4f} u + {gCal:.4f} phi  "
      f"(R^2 = {calR2:.4f}, rms = {np.std(calRes):.4f})")
print(f"state radius R = {rMed:.3f}, logged axis length |(c,s)| = {axisMag:.3f}")
bEff = bCal * rMed * axisMag  # cos(D) coefficient as this log was run
print(f"cos(D) coefficient / constant: this log {bEff / a0Cal:.3f}, "
      f"with unit axis {bCal * rMed / a0Cal:.3f}  (ideal 1)")
phiLoRec = (a0Cal - bCal * rMed) / abs(gCal)
phiHiRec = (a0Cal + bCal * rMed) / abs(gCal)
print(f"Recommended phi sweep (unit axis): [{phiLoRec:.4f}, {phiHiRec:.4f}] of 5 V "
      f"-> V_PHI_LO = {5 * phiLoRec:.3f} V, V_PHI_HI = {5 * phiHiRec:.3f} V")
phiLoLog, phiHiLog = float(np.min(ph[calIdx])), float(np.max(ph[calIdx]))
print(f"This log swept phi over [{phiLoLog:.4f}, {phiHiLog:.4f}] "
      f"({len(np.unique(np.round(ph[calIdx], 4)))} values)")

# Predicted P(up) for this log's settings: continuous sweep and the actual phi values
def predictUp(Dgrid, phiVals):
    thr = (a0Cal + bEff * np.cos(Dgrid))[:, None] + gCal * phiVals[None, :]
    return np.mean(thr > 0, axis = 1)

phiCont = np.linspace(phiLoLog, phiHiLog, 4001)
phiDisc = np.unique(ph[calIdx])
predCont = predictUp(thetaDom, phiCont)
predDiscBins = predictUp(binThets, phiDisc)
predContBins = predictUp(binThets, phiCont)
idealBins = np.cos(0.5 * binThets)**2
rmsData = np.sqrt(np.nanmean((probs - idealBins)**2))
rmsCal = np.sqrt(np.mean((predContBins - idealBins)**2))
rmsDisc = np.sqrt(np.nanmean((probs - predDiscBins)**2))
print(f"RMS deviation from Born rule: data {rmsData:.3f}; "
      f"calibration alone (continuous phi) {rmsCal:.3f}; "
      f"data vs discrete-phi prediction {rmsDisc:.3f}")

calFig, calAx = plt.subplots()
calAx.plot(thetaDom, ideal, label = 'Born rule')
calAx.plot(thetaDom, predCont, 'k--', lw = 1,
           label = 'predicted, continuous $\\phi$ (this calibration)')
calAx.plot(binThets, predDiscBins, 's', mfc = 'none', color = 'gray',
           label = f'predicted, {len(phiDisc)} $\\phi$ values')
calAx.plot(binThets[valid], probs[valid], '.', color = 'red', markersize = 14,
           label = 'empirical')
calAx.set_xlabel("Relative angle (rad)")
calAx.set_ylabel("Up probability")
calAx.set_title('Empirical vs calibration-predicted up probability')
calAx.legend(fontsize = 8)
calFig.tight_layout()
calFig.savefig('bornCalibration.pdf')
calFig.show()


'''
ODE FIDELITY: does (x, y) collapse to sign(theta) * (c, s)?

Measurement force is the on/off column (of), not (c, s) leaving (0, 0).
Within each cycle, (c, s) stays fixed at the axis while of is True for the
first half (force on) and False for the second; (c, s) = (0, 0) only appears
at the start of the log before the first measurement.

Predictor: sign of logged theta (corrected conversion) at the first force-on
sample after c/s have settled (thetaLag samples in for pre-v2 logs), not the
mean over the collapse window.
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

    # Theta at the start of the collapse, once c/s have settled (see thetaLag)
    th0 = th[min(startIdx + thetaLag, forceIdx)]
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
print(f"Predictor: sign of theta {thetaLag} samples after force-on start")
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

Logged x, y are ADC fractions of 5 V with observed |r| ≈ 0.87 (the 14/19
bias in the 1-r^2 circuit), so (1 - r_ADC^2) does not vanish on the physical
circle. Rescale
  (x, y)_nat = (x, y)_ADC / R_obs
with R_obs = median |r|_ADC, then fit in natural units:
  a ≈ k q (1 - r_nat^2) - d v + g theta n_eff
Theta and c, s keep their own voltage scales; g absorbs that mismatch.

v is not logged: estimate v and a = dv/dt from finite differences in t.
Firmware: when of is False, measureOff() writes PWM 128 (≈ +0.0039 of 5 V),
but the log still stores the software axis. Use n_eff = (c, s) when on and
(C_OFF, C_OFF) when off.

CAVEAT: in the circuit the velocity stage (1 kΩ || 10 nF) has a 10 µs time
constant, so v follows its inputs instantly and the position integrator
(1 MΩ, 10 nF) sets a 10 ms time scale. At ~60 ms sampling these
acceleration-based fits cannot resolve the dynamics; they are kept for
reference. See the FIRST-ORDER sections at the end for the fits that match
the circuit.
'''

# Observed circle radius in logged (ADC fraction) units
rAdc = np.sqrt(x**2 + y**2)
Robs = float(np.median(rAdc))
xN, yN = x / Robs, y / Robs
# c, s are commanded drive voltages (fraction of 5V); the 14/19 radius
# artifact lives in the 1-r^2 circuit, so c, s are NOT rescaled by R_obs.

# Effective measurement axis on the wires: commanded axis when on,
# the PWM-128 offset when off
cEff = np.where(of, c, C_OFF)
sEff = np.where(of, s, C_OFF)

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
ONE-STEP PREDICTION + REDUCED CHI-SQUARED (FULL SECOND-ORDER STATE)

The ODEs are second order, so the state is (x, y, v_x, v_y). From the measured
state at t_i, predict the full state at t_{i+1} and score all four residuals.

  chi^2 = sum_i sum_{alpha in {x,y,vx,vy}} (r_{alpha,i} / sigma_alpha)^2
  nu    = 4 N_steps - 3   (subtract fitted k, d, g)

One-step Euler map (linear in k, d, g):
  a = k q (1-r^2) - d v + g theta n_eff
  v_{i+1} = v_i + dt a
  q_{i+1} = q_i + dt v_i + dt^2 a

Velocities are causal backward differences of the smoothed positions (no look-ahead).
Sigmas are force-off one-step residuals after fitting (k, d) with n_eff = 0.
'''

print("\n--- One-step ODE check, full state (x,y,vx,vy) ---")

# Causal velocity from smoothed natural-unit positions
vxC = np.full(len(t), np.nan)
vyC = np.full(len(t), np.nan)
vxC[1:] = np.diff(xSoft) / np.diff(t)
vyC[1:] = np.diff(ySoft) / np.diff(t)

# Pairs need a defined v at i (so i >= 1) and valid smoothing support
offPairs = []
onPairs = []
for i in range(len(t) - 1):
    if not (validSoft[i] and validSoft[i + 1]):
        continue
    if not np.isfinite(vxC[i]):
        continue
    if (not of[i]) and (not of[i + 1]):
        offPairs.append(i)
    elif of[i] and of[i + 1]:
        onPairs.append(i)
offPairs = np.array(offPairs, dtype = int)
onPairs = np.array(onPairs, dtype = int)

# Channel order within each time step: x, y, vx, vy
STATE_LABELS = ('x', 'y', 'vx', 'vy')

def oneStepDesignFull(pairs):
    """
    Four rows per time step (x, y, vx, vy).
    Position target: q_{i+1} - q_i - dt v_i   ~ dt^2 * a
    Velocity target: v_{i+1} - v_i           ~ dt  * a
    with a columns [q(1-r^2), -v, theta n].
    """
    rows = []
    targets = []
    for i in pairs:
        dtStep = t[i + 1] - t[i]
        dt2 = dtStep * dtStep
        r2i = xSoft[i]**2 + ySoft[i]**2
        radX = xSoft[i] * (1.0 - r2i)
        radY = ySoft[i] * (1.0 - r2i)
        # x position
        rows.append((dt2 * radX, dt2 * (-vxC[i]), dt2 * thSoft[i] * cEff[i]))
        targets.append(xSoft[i + 1] - xSoft[i] - dtStep * vxC[i])
        # y position
        rows.append((dt2 * radY, dt2 * (-vyC[i]), dt2 * thSoft[i] * sEff[i]))
        targets.append(ySoft[i + 1] - ySoft[i] - dtStep * vyC[i])
        # vx
        rows.append((dtStep * radX, dtStep * (-vxC[i]), dtStep * thSoft[i] * cEff[i]))
        targets.append(vxC[i + 1] - vxC[i])
        # vy
        rows.append((dtStep * radY, dtStep * (-vyC[i]), dtStep * thSoft[i] * sEff[i]))
        targets.append(vyC[i + 1] - vyC[i])
    return np.asarray(rows, dtype = float), np.asarray(targets, dtype = float)

def splitChannels(resid):
    """Split flat residual vector into dict of per-channel arrays."""
    return {lab: resid[i::4] for i, lab in enumerate(STATE_LABELS)}

def distMoments(arr, name):
    arr = np.asarray(arr, dtype = float)
    mu = float(np.mean(arr))
    std = float(np.std(arr))
    centered = arr - mu
    skew = float(np.mean(centered**3) / std**3) if std > 0 else float('nan')
    kurt = float(np.mean(centered**4) / std**4) if std > 0 else float('nan')
    mad = float(np.median(np.abs(arr - np.median(arr))))
    sigmaMad = 1.4826 * mad
    z = np.abs(centered) / std if std > 0 else np.zeros_like(arr)
    print(f"{name}: n={len(arr)}, mean={mu:.4g}, std={std:.4g}, "
          f"skew={skew:.3f}, kurtosis={kurt:.3f} (Gauss=3, excess={kurt - 3:.3f})")
    print(f"  MAD={mad:.4g}, sigma_MAD={sigmaMad:.4g}, "
          f"std/sigma_MAD={std / sigmaMad if sigmaMad > 0 else float('nan'):.2f}")
    print(f"  frac |z|>2,3,5: {np.mean(z > 2):.3f}, {np.mean(z > 3):.3f}, "
          f"{np.mean(z > 5):.3f}  (Gauss ~0.046, 0.0027, 6e-7)")
    return dict(mu = mu, std = std, skew = skew, kurt = kurt, sigmaMad = sigmaMad)

# --- Force-off: fit k, d; sigma per channel from residuals ---
AoffStep, boffStep = oneStepDesignFull(offPairs)
AoffKD = AoffStep[:, :2]
coeffOff, _, _, _ = np.linalg.lstsq(AoffKD, boffStep, rcond = None)
kOff, dOff = coeffOff
residOff = boffStep - AoffKD @ coeffOff
chOff = splitChannels(residOff)
sigmas = {lab: float(np.sqrt(np.mean(chOff[lab]**2))) for lab in STATE_LABELS}

nOffObs = len(boffStep)
nOffParams = 2
nuOffSigma = nOffObs - nOffParams
chi2OffCal = sum(np.sum((chOff[lab] / sigmas[lab])**2) for lab in STATE_LABELS)

adcStepNat = logStep / Robs
print(f"Force-off one-step pairs: {len(offPairs)}")
print(f"Force-off ODE fit (k, d only): k={kOff:.3g}, d={dOff:.3g}")
print(f"Tight sigmas (q and v) after force-off k,d fit:")
for lab in STATE_LABELS:
    print(f"  sigma_{lab} = {sigmas[lab]:.4g}")
print(f"Logged position step (natural units): {adcStepNat:.4g}")
print(f"Force-off chi^2_red (calibration): {chi2OffCal / nuOffSigma:.3f}  "
      f"(nu={nuOffSigma})")

print("Force-off residual distributions (sigma samples):")
momsOff = {lab: distMoments(chOff[lab], f"  d{lab}") for lab in STATE_LABELS}

sigFig, sigAxes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))
zg = np.linspace(-6, 6, 300)
gauss = np.exp(-0.5 * zg**2) / np.sqrt(2 * np.pi)
for axS, lab in zip(sigAxes.ravel(), STATE_LABELS):
    z = chOff[lab] / sigmas[lab]
    axS.hist(z, bins = 80, density = True, alpha = 0.7, color = 'tab:blue')
    axS.plot(zg, gauss, 'k--', lw = 1)
    axS.set_xlim(-6, 6)
    axS.set_title(lab)
    axS.set_xlabel(r'residual / $\sigma$')
sigFig.suptitle('Force-off full-state residuals (sigma samples)')
sigFig.tight_layout()
sigFig.savefig('odeSigmaDistribution.pdf')
sigFig.show()

# --- Force-on: weighted fit of k, d, g on full state ---
AonStep, bonStep = oneStepDesignFull(onPairs)
w = np.empty(len(bonStep))
for i, lab in enumerate(STATE_LABELS):
    w[i::4] = 1.0 / sigmas[lab]
coeffStep, _, _, _ = np.linalg.lstsq(AonStep * w[:, None], bonStep * w, rcond = None)
kStep, dStep, gStep = coeffStep
residOn = bonStep - AonStep @ coeffStep
chOn = splitChannels(residOn)

nParams = 3
nStepsOn = len(onPairs)
nObsOn = 4 * nStepsOn
nuOn = nObsOn - nParams
chi2On = sum(np.sum((chOn[lab] / sigmas[lab])**2) for lab in STATE_LABELS)
chi2RedOn = chi2On / nuOn if nuOn > 0 else float('nan')

# Position-only and velocity-only chi^2 for diagnosis
chi2Pos = (np.sum((chOn['x'] / sigmas['x'])**2) + np.sum((chOn['y'] / sigmas['y'])**2))
chi2Vel = (np.sum((chOn['vx'] / sigmas['vx'])**2) + np.sum((chOn['vy'] / sigmas['vy'])**2))
nuPos = 2 * nStepsOn  # not subtracting params again for the split view
nuVel = 2 * nStepsOn

print(f"One-step fit on force-on pairs (N_steps={nStepsOn}):")
print(f"  k={kStep:.3g}, d={dStep:.3g}, g={gStep:.3g}")
print(f"  chi^2 = {chi2On:.3g},  nu = 4*{nStepsOn} - {nParams} = {nuOn}")
print(f"  reduced chi^2 (full state) = {chi2RedOn:.3f}")
print(f"  position-only chi^2 / 2N = {chi2Pos / nuPos:.3f}")
print(f"  velocity-only chi^2 / 2N = {chi2Vel / nuVel:.3f}")
print("  RMS residuals vs sigma:")
for lab in STATE_LABELS:
    rms = float(np.sqrt(np.mean(chOn[lab]**2)))
    print(f"    {lab}: rms={rms:.4g}, sigma={sigmas[lab]:.4g}, ratio={rms / sigmas[lab]:.2f}")

print("Force-on residual distributions:")
for lab in STATE_LABELS:
    distMoments(chOn[lab], f"  d{lab}")

chi2PerStep = sum((chOn[lab] / sigmas[lab])**2 for lab in STATE_LABELS)
print(f"Per-step chi^2: mean={np.mean(chi2PerStep):.2f}, median={np.median(chi2PerStep):.2f}, "
      f"p95={np.percentile(chi2PerStep, 95):.2f}, max={np.max(chi2PerStep):.2f}")
top5 = chi2PerStep >= np.percentile(chi2PerStep, 95)
print(f"Fraction of total chi^2 from top 5% steps: "
      f"{np.sum(chi2PerStep[top5]) / np.sum(chi2PerStep):.2f}")

# Heteroscedasticity in q and v
stepQOff = np.hypot(xSoft[offPairs + 1] - xSoft[offPairs],
                    ySoft[offPairs + 1] - ySoft[offPairs])
stepQOn = np.hypot(xSoft[onPairs + 1] - xSoft[onPairs],
                   ySoft[onPairs + 1] - ySoft[onPairs])
stepVOff = np.hypot(vxC[offPairs + 1] - vxC[offPairs],
                    vyC[offPairs + 1] - vyC[offPairs])
stepVOn = np.hypot(vxC[onPairs + 1] - vxC[onPairs],
                   vyC[onPairs + 1] - vyC[onPairs])
print(f"RMS |Delta q| per step: off={np.sqrt(np.mean(stepQOff**2)):.4g}, "
      f"on={np.sqrt(np.mean(stepQOn**2)):.4g}, "
      f"ratio={np.sqrt(np.mean(stepQOn**2) / np.mean(stepQOff**2)):.2f}")
print(f"RMS |Delta v| per step: off={np.sqrt(np.mean(stepVOff**2)):.4g}, "
      f"on={np.sqrt(np.mean(stepVOn**2)):.4g}, "
      f"ratio={np.sqrt(np.mean(stepVOn**2) / np.mean(stepVOff**2)):.2f}")

# Explanatory power: ballistic/coasting baseline (a=0) vs ODE
# coasting targets are exactly bonStep; ODE residual is residOn
baseVar = float(np.mean(bonStep**2))
odeVar = float(np.mean(residOn**2))
print(f"Fraction of coasting one-step variance explained by ODE (all channels): "
      f"{1.0 - odeVar / baseVar:.3f}")

keep = ~top5
nuTrim = int(np.sum(keep)) * 4 - nParams
print(f"Chi^2_red after dropping worst 5% steps: "
      f"{np.sum(chi2PerStep[keep]) / nuTrim:.3f} (nu={nuTrim})")

# Residual histograms by channel
chiFig, chiAxes = plt.subplots(nrows = 2, ncols = 2, figsize = (8, 6))
for axC, lab in zip(chiAxes.ravel(), STATE_LABELS):
    axC.hist(chOn[lab] / sigmas[lab], bins = 60, density = True, alpha = 0.7)
    axC.plot(zg, gauss, 'k--', lw = 1)
    axC.set_xlim(-6, 6)
    axC.set_title(lab)
    axC.set_xlabel(r'residual / $\sigma$')
chiFig.suptitle(r'Force-on full-state residuals ($\chi^2_\nu = %0.2f$)' % chi2RedOn)
chiFig.tight_layout()
chiFig.savefig('odeOneStepResiduals.pdf')
chiFig.show()

# Scatter: predicted vs measured increment for position and velocity (x-channels)
stepFig, stepAxes = plt.subplots(nrows = 1, ncols = 2, figsize = (8, 4))
predOn = AonStep @ coeffStep
for axP, idx, title in (
    (stepAxes[0], 0, r'position $x$: $q_{i+1}-q_i-\Delta t v_i$'),
    (stepAxes[1], 2, r'velocity $v_x$: $v_{i+1}-v_i$'),
):
    meas = bonStep[idx::4]
    pred = predOn[idx::4]
    nShowStep = min(3000, len(meas))
    showStep = rng.choice(len(meas), size = nShowStep, replace = False)
    axP.scatter(meas[showStep], pred[showStep], s = 6, alpha = 0.25)
    limStep = np.percentile(np.abs(meas), 99)
    axP.plot([-limStep, limStep], [-limStep, limStep], 'k--', lw = 1)
    axP.set_xlabel('measured')
    axP.set_ylabel('ODE')
    axP.set_title(title)
    axP.set_aspect('equal', adjustable = 'box')
    axP.set_xlim(-limStep, limStep)
    axP.set_ylim(-limStep, limStep)
stepFig.suptitle('One-step full-state force fit')
stepFig.tight_layout()
stepFig.savefig('odeOneStepScatter.pdf')
stepFig.show()


'''
FORCE LIBRARY IDENTIFICATION (what is the analog computer actually computing?)

Estimate a ≈ Δv/Δt from the causal one-step state, then regress onto a library of
plausible analog-computer terms built from (x, y, v, c_eff, s_eff, theta, r^2, u).
Mirrored libraries for a_x and a_y share one coefficient vector.

Designed ODE terms are included explicitly:
  x(1-r^2), -v_x, theta*c   (and y-mirrors)
plus extras that multipliers/summers could produce.
'''

print("\n--- Force library identification ---")

LIB_NAMES = [
    '1',
    'x', 'y',                 # with x*r^2 spans designed x*(1-r^2)
    'vx', 'vy',
    'c', 's',
    'theta',
    'theta*c', 'theta*s',     # designed measurement force
    'theta*x', 'theta*y',
    'x*r^2', 'y*r^2',
    'u',                      # u = x c + y s
    'u*c', 'u*s',             # appears if force ~ (1+u)*axis
    'u*x', 'u*y',
    'x*c', 'x*s', 'y*c', 'y*s',
    'theta*u',
    'theta*u*c', 'theta*u*s',
]

def featuresDict(qx, qy, vx, vy, cx, sy, thv):
    r2v = qx * qx + qy * qy
    u = qx * cx + qy * sy
    return {
        '1': 1.0,
        'x': qx, 'y': qy,
        'vx': vx, 'vy': vy,
        'c': cx, 's': sy,
        'theta': thv,
        'theta*c': thv * cx, 'theta*s': thv * sy,
        'theta*x': thv * qx, 'theta*y': thv * qy,
        'x*r^2': qx * r2v, 'y*r^2': qy * r2v,
        'u': u,
        'u*c': u * cx, 'u*s': u * sy,
        'u*x': u * qx, 'u*y': u * qy,
        'x*c': qx * cx, 'x*s': qx * sy, 'y*c': qy * cx, 'y*s': qy * sy,
        'theta*u': thv * u,
        'theta*u*c': thv * u * cx, 'theta*u*s': thv * u * sy,
    }

# Shared coeffs: a_y uses mirrored features under x<->y, c<->s, vx<->vy
MIRROR_MAP = {
    '1': '1',
    'x': 'y', 'y': 'x',
    'vx': 'vy', 'vy': 'vx',
    'c': 's', 's': 'c',
    'theta': 'theta',
    'theta*c': 'theta*s', 'theta*s': 'theta*c',
    'theta*x': 'theta*y', 'theta*y': 'theta*x',
    'x*r^2': 'y*r^2', 'y*r^2': 'x*r^2',
    'u': 'u',
    'u*c': 'u*s', 'u*s': 'u*c',
    'u*x': 'u*y', 'u*y': 'u*x',
    'x*c': 'y*s', 'y*s': 'x*c',
    'x*s': 'y*c', 'y*c': 'x*s',
    'theta*u': 'theta*u',
    'theta*u*c': 'theta*u*s', 'theta*u*s': 'theta*u*c',
}

def buildForceSystem(pairs):
    rows = []
    targets = []
    for i in pairs:
        dtStep = t[i + 1] - t[i]
        if dtStep <= 0 or not np.isfinite(vxC[i + 1]):
            continue
        f = featuresDict(
            xSoft[i], ySoft[i], vxC[i], vyC[i],
            cEff[i], sEff[i], thSoft[i],
        )
        rows.append([f[n] for n in LIB_NAMES])
        targets.append((vxC[i + 1] - vxC[i]) / dtStep)
        rows.append([f[MIRROR_MAP[n]] for n in LIB_NAMES])
        targets.append((vyC[i + 1] - vyC[i]) / dtStep)
    return np.asarray(rows, dtype = float), np.asarray(targets, dtype = float)

def sindyFit(A, b, names, nThresh = 8, threshFrac = 0.05):
    """Least squares with iterative hard thresholding (simple SINDy)."""
    coeff = np.linalg.lstsq(A, b, rcond = None)[0]
    for _ in range(nThresh):
        scale = np.max(np.abs(coeff)) if np.max(np.abs(coeff)) > 0 else 1.0
        keep = np.abs(coeff) >= threshFrac * scale
        if np.sum(keep) == 0:
            keep[np.argmax(np.abs(coeff))] = True
        coeff2 = np.zeros_like(coeff)
        coeff2[keep] = np.linalg.lstsq(A[:, keep], b, rcond = None)[0]
        coeff = coeff2
    return coeff

def reportForceFit(pairs, label, doSindy = True, dropAxis = False):
    A, b = buildForceSystem(pairs)
    # Drop near-collinear empty columns (e.g. c,s during force-off)
    colStd = np.std(A, axis = 0)
    active = colStd > 1e-9
    if dropAxis:
        # Force-off: the axis is only the tiny constant PWM offset, so terms
        # built from c, s, u are either duplicates or scaled-up noise fits.
        for j, n in enumerate(LIB_NAMES):
            if set(n.split('*')) & {'c', 's', 'u'}:
                active[j] = False
    # Also drop columns that are (near) linear combinations of earlier ones.
    # During force-off, c = s = C_OFF is constant, so e.g. theta*c is just
    # C_OFF*theta and u*x duplicates x^2 + xy terms; keeping them makes the
    # least-squares coefficients meaningless.
    kept = []
    for j in np.where(active)[0]:
        trial = kept + [j]
        Aj = A[:, trial] / np.linalg.norm(A[:, trial], axis = 0)
        sv = np.linalg.svd(Aj, compute_uv = False)
        if sv[-1] / sv[0] > 1e-6:
            kept.append(j)
        else:
            active[j] = False
    Aact = A[:, active]
    namesAct = [n for n, a in zip(LIB_NAMES, active) if a]

    coeffFull = np.zeros(len(LIB_NAMES))
    coeffFull[active] = np.linalg.lstsq(Aact, b, rcond = None)[0]
    predFull = A @ coeffFull
    r2Full = rSquared(b, predFull)
    rms = float(np.sqrt(np.mean((b - predFull)**2)))

    print(f"\n{label}: N_eq={len(b)} ({len(pairs)} steps x 2)")
    print(f"  dense library R^2={r2Full:.3f}, residual RMS={rms:.4g}")
    order = np.argsort(-np.abs(coeffFull))
    print("  Top |coeff| terms (dense):")
    for j in order[:12]:
        if abs(coeffFull[j]) < 1e-9:
            continue
        print(f"    {LIB_NAMES[j]:16s}  {coeffFull[j]:+.4g}")

    # Designed ODE subspace: a = k x(1-r^2) - d v + g theta*c
    # in this library: k*x + (-k)*x*r^2 + (-d)*vx + g*theta*c
    designed = ['x', 'x*r^2', 'vx', 'theta*c']
    dIdx = [LIB_NAMES.index(n) for n in designed if active[LIB_NAMES.index(n)]]
    if len(dIdx) >= 3:
        Ad = A[:, dIdx]
        cd = np.linalg.lstsq(Ad, b, rcond = None)[0]
        predD = Ad @ cd
        r2D = rSquared(b, predD)
        print(f"  Designed ODE subspace {designed}: R^2={r2D:.3f}")
        for n, c in zip([LIB_NAMES[i] for i in dIdx], cd):
            print(f"    {n:16s}  {c:+.4g}")
        nameToC = {LIB_NAMES[i]: c for i, c in zip(dIdx, cd)}
        if 'x' in nameToC and 'x*r^2' in nameToC:
            print(f"    implied k from 'x': {nameToC['x']:+.4g}, "
                  f"from '-coeff(x*r^2)': {-nameToC['x*r^2']:+.4g}")
        if 'vx' in nameToC:
            print(f"    implied d = -coeff(vx): {-nameToC['vx']:+.4g}")
        if 'theta*c' in nameToC:
            print(f"    implied g = coeff(theta*c): {nameToC['theta*c']:+.4g}")

    if doSindy:
        coeffSp = np.zeros(len(LIB_NAMES))
        coeffSp[active] = sindyFit(Aact, b, namesAct)
        predSp = A @ coeffSp
        r2Sp = rSquared(b, predSp)
        nnz = int(np.sum(np.abs(coeffSp) > 0))
        print(f"  SINDy-sparse R^2={r2Sp:.3f}, {nnz} active terms:")
        orderSp = np.argsort(-np.abs(coeffSp))
        for j in orderSp:
            if abs(coeffSp[j]) < 1e-12:
                continue
            print(f"    {LIB_NAMES[j]:16s}  {coeffSp[j]:+.4g}")
        return coeffFull, coeffSp, r2Full, r2Sp

    return coeffFull, None, r2Full, None

coeffOnDense, coeffOnSp, r2OnDense, r2OnSp = reportForceFit(
    onPairs, "Force-on library fit", doSindy = True,
)
coeffOffDense, coeffOffSp, r2OffDense, r2OffSp = reportForceFit(
    offPairs, "Force-off library fit", doSindy = True, dropAxis = True,
)

# Coefficient comparison bar chart (force-on dense, top terms)
libFig, libAx = plt.subplots(figsize = (9, 4))
orderPlot = np.argsort(-np.abs(coeffOnDense))[:15]
namesPlot = [LIB_NAMES[j] for j in orderPlot]
valsPlot = coeffOnDense[orderPlot]
colors = ['tab:red' if n in ('theta*c', 'x', 'x*r^2', 'vx') else 'tab:blue'
          for n in namesPlot]
libAx.barh(range(len(namesPlot))[::-1], valsPlot[::-1], color = colors[::-1])
libAx.set_yticks(range(len(namesPlot)))
libAx.set_yticklabels(namesPlot[::-1])
libAx.axvline(0.0, color = 'k', lw = 0.8)
libAx.set_xlabel('Coefficient in a ≈ sum c_k f_k')
libAx.set_title(f'Force-on library terms (dense fit, $R^2={r2OnDense:.2f}$; '
                f'red = designed ODE)')
libFig.tight_layout()
libFig.savefig('odeForceLibrary.pdf')
libFig.show()


'''
HYPOTHESIS: measurement force is g (1 + u - phi) n, not g theta_logged n

u = x c + y s (natural units). Compare nested one-step acceleration models
on force-on data, and check whether logged theta matches 1 + u - phi.
'''

print("\n--- Hypothesis: force ~ (1 + u - phi) * n ---")

# phi is commanded/stored per sample; use soft-averaged only if we want — keep raw ph
# aligned with causal state at index i
phiSoft = movingAverage(ph, softWin)  # mild smooth for consistency with other channels

def hypTargets(pairs):
    """Measured a_x, a_y from causal Delta v / dt."""
    axm, aym = [], []
    for i in pairs:
        dtStep = t[i + 1] - t[i]
        if dtStep <= 0 or not np.isfinite(vxC[i + 1]):
            continue
        axm.append((vxC[i + 1] - vxC[i]) / dtStep)
        aym.append((vyC[i + 1] - vyC[i]) / dtStep)
    return np.asarray(axm), np.asarray(aym)

def hypDesign(pairs, thetaMode, phiScale = 1.0, phiOffset = 0.0):
    """
    Columns for stacked [a_x; a_y]:
      [x(1-r^2), -vx,  f * c ] and mirrored for y
    where f depends on thetaMode.
    """
    rows = []
    targets = []
    extras = []  # store f used, u, phi, theta for diagnostics
    for i in pairs:
        dtStep = t[i + 1] - t[i]
        if dtStep <= 0 or not np.isfinite(vxC[i + 1]):
            continue
        qx, qy = xSoft[i], ySoft[i]
        vx, vy = vxC[i], vyC[i]
        cx, sy = cEff[i], sEff[i]
        r2v = qx * qx + qy * qy
        u = qx * cx + qy * sy
        phiV = phiSoft[i] * phiScale + phiOffset
        thLog = thSoft[i]

        if thetaMode == 'logged':
            f = thLog
        elif thetaMode == '1+u-phi':
            f = 1.0 + u - phiV
        elif thetaMode == 'u-phi':
            f = u - phiV
        elif thetaMode == '1+u':
            f = 1.0 + u
        else:
            raise ValueError(thetaMode)

        # a_x row
        rows.append((qx * (1.0 - r2v), -vx, f * cx))
        targets.append((vxC[i + 1] - vxC[i]) / dtStep)
        # a_y row
        rows.append((qy * (1.0 - r2v), -vy, f * sy))
        targets.append((vyC[i + 1] - vyC[i]) / dtStep)
        extras.append((f, u, phiV, thLog, cx, sy))

    return (np.asarray(rows, dtype = float),
            np.asarray(targets, dtype = float),
            extras)

def fitHyp(pairs, thetaMode, label, phiScale = 1.0, phiOffset = 0.0):
    A, b, extras = hypDesign(pairs, thetaMode, phiScale, phiOffset)
    coeff, _, _, _ = np.linalg.lstsq(A, b, rcond = None)
    pred = A @ coeff
    r2 = rSquared(b, pred)
    kH, dH, gH = coeff
    print(f"{label}:")
    print(f"  k={kH:.4g}, d={dH:.4g}, g={gH:.4g},  R^2={r2:.3f}, "
          f"RMS={np.sqrt(np.mean((b - pred)**2)):.4g}")
    return dict(coeff = coeff, r2 = r2, pred = pred, meas = b, extras = extras, A = A)

# Compare core hypotheses on force-on
print("Force-on nested models a = k q(1-r^2) - d v + g f n:")
fitLog = fitHyp(onPairs, 'logged', "  f = theta_logged")
fitHyp1 = fitHyp(onPairs, '1+u-phi', "  f = 1 + u - phi_logged")
fitHyp2 = fitHyp(onPairs, 'u-phi', "  f = u - phi_logged")
fitHyp3 = fitHyp(onPairs, '1+u', "  f = 1 + u  (ignore phi)")

# Free linear f = b0 + b1 u + b2 phi  (does data want 1, +u, -phi?)
rowsF, targsF = [], []
for i in onPairs:
    dtStep = t[i + 1] - t[i]
    if dtStep <= 0 or not np.isfinite(vxC[i + 1]):
        continue
    qx, qy = xSoft[i], ySoft[i]
    vx, vy = vxC[i], vyC[i]
    cx, sy = cEff[i], sEff[i]
    r2v = qx * qx + qy * qy
    u = qx * cx + qy * sy
    phiV = phiSoft[i]
    # a_x: [q(1-r2), -v, 1*n, u*n, phi*n]
    rowsF.append((qx * (1.0 - r2v), -vx, cx, u * cx, phiV * cx))
    targsF.append((vxC[i + 1] - vxC[i]) / dtStep)
    rowsF.append((qy * (1.0 - r2v), -vy, sy, u * sy, phiV * sy))
    targsF.append((vyC[i + 1] - vyC[i]) / dtStep)
AF = np.asarray(rowsF)
bF = np.asarray(targsF)
cF = np.linalg.lstsq(AF, bF, rcond = None)[0]
predF = AF @ cF
r2F = rSquared(bF, predF)
print("  f = b0 + b1 u + b2 phi  (free), via g*f*n expanded:")
print(f"    k={cF[0]:.4g}, d={cF[1]:.4g}, g*b0={cF[2]:.4g}, g*b1={cF[3]:.4g}, "
      f"g*b2={cF[4]:.4g},  R^2={r2F:.3f}")
if abs(cF[2]) > 1e-12:
    print(f"    implied (b0, b1, b2) ~ (1, {cF[3]/cF[2]:.3f}, {cF[4]/cF[2]:.3f}) "
          f"if normalized to b0=1;  hypothesis (1,+1,-1)")

# Try phi rescaled to roughly [0, 2] as in the math ODE
phiMaxLog = float(np.max(ph))
phiMinLog = float(np.min(ph[ph > 0])) if np.any(ph > 0) else 0.0
print(f"\nLogged phi range: [{np.min(ph):.4g}, {np.max(ph):.4g}]")
for scaleName, scale in [
    ('phi as logged', 1.0),
    ('phi * (2/phiMax)', 2.0 / phiMaxLog if phiMaxLog > 0 else 1.0),
    ('phi * 5 (volts)', 5.0),
    ('phi * (2/0.36)', 2.0 / 0.36),
]:
    fitHyp(onPairs, '1+u-phi', f"  f=1+u-phi with {scaleName}", phiScale = scale)

# Does theta_logged track 1+u-phi?
uOn = []
phiOn = []
thOn = []
hypOn = []
for i in onPairs:
    cx, sy = cEff[i], sEff[i]
    u = xSoft[i] * cx + ySoft[i] * sy
    uOn.append(u)
    phiOn.append(phiSoft[i])
    thOn.append(thSoft[i])
    hypOn.append(1.0 + u - phiSoft[i])
uOn, phiOn, thOn, hypOn = map(np.asarray, (uOn, phiOn, thOn, hypOn))
corr = np.corrcoef(thOn, hypOn)[0, 1]
print(f"\ntheta_logged vs 1+u-phi_logged: corr={corr:.3f}, "
      f"RMSE={np.sqrt(np.mean((thOn - hypOn)**2)):.3f}")
print(f"  mean theta_logged={np.mean(thOn):.3f}, mean(1+u-phi)={np.mean(hypOn):.3f}")
print(f"  best affine theta ≈ a + b(1+u-phi): ", end = '')
Ba = np.column_stack((np.ones_like(hypOn), hypOn))
ab = np.linalg.lstsq(Ba, thOn, rcond = None)[0]
print(f"a={ab[0]:.3f}, b={ab[1]:.3f}, R^2={rSquared(thOn, Ba @ ab):.3f}")

# Scatter theta vs hypothesis
thFig, thAx = plt.subplots()
thAx.scatter(hypOn, thOn, s = 6, alpha = 0.2, label = 'force-on samples')
lims = [min(hypOn.min(), thOn.min()), max(hypOn.max(), thOn.max())]
thAx.plot(lims, lims, 'k--', lw = 1, label = 'equality')
thAx.plot(lims, ab[0] + ab[1] * np.asarray(lims), 'r-', lw = 1, label = 'affine fit')
thAx.set_xlabel(r'$1 + u - \phi_{\mathrm{logged}}$')
thAx.set_ylabel(r'$\theta_{\mathrm{logged}}$')
thAx.set_title(r'Does logged $\theta$ match $1+u-\phi$?')
thAx.legend()
thFig.tight_layout()
thFig.savefig('odeThetaHypothesis.pdf')
thFig.show()

# Model comparison bar chart
hypFig, hypAx = plt.subplots()
namesH = ['theta_logged', '1+u-phi', 'u-phi', '1+u', 'free b0+b1 u+b2 phi']
r2sH = [fitLog['r2'], fitHyp1['r2'], fitHyp2['r2'], fitHyp3['r2'], r2F]
hypAx.barh(namesH[::-1], r2sH[::-1], color = 'tab:purple')
hypAx.set_xlabel(r'$R^2$ of acceleration model')
hypAx.set_title(r'Force model: $a=kq(1-r^2)-dv+g f\,n$')
hypFig.tight_layout()
hypFig.savefig('odeForceHypothesisCompare.pdf')
hypFig.show()


'''
FIRST-ORDER TANGENTIAL MODEL (matches the circuit time scales)

With v slaved to its inputs, the circuit is first order:
  dq/dt = (1/tau) * [radial(q) + (theta c, theta s) + F0]
tau = 1 MΩ * 10 nF = 10 ms. At 60 ms sampling the radial part is always at
equilibrium, so score the TANGENTIAL motion, where the radial force drops out:
  r dpsi/dt = g theta (n_eff . t_hat) + F0 . t_hat,
  t_hat = (-sin psi, cos psi)
All quantities are in logged units (fraction of 5 V). Predicted from the
schematic: AD633 W = V_theta V_c / 10 V, unity-gain velocity stage, so
  g_pred = 5 V / (10 V * tau) = 50 /s,
and a velocity-node offset V_off gives F0 = V_off / (5 V * tau).

One fit covers both windows: force-on uses the commanded axis, force-off uses
the PWM-128 offset (C_OFF, C_OFF). F0 is the drift from op-amp/multiplier
offsets that is present whether or not a measurement is running.
'''

print("\n--- First-order tangential model (drift + measurement force) ---")
TAU_PRED = 10e-3
G_PRED = 5.0 / (10.0 * TAU_PRED)

psiL = np.arctan2(y, x)
rL = np.hypot(x, y)

def tangentialRows(pairs):
    """Midpoint-rule rows [theta n.t, t_x, t_y] and target r dpsi/dt."""
    i = np.asarray(pairs, dtype = int)
    j = i + 1
    dtP = t[j] - t[i]
    dpsi = wrapAngle(psiL[j] - psiL[i])
    psiM = psiL[i] + 0.5 * dpsi
    tx, ty = -np.sin(psiM), np.cos(psiM)
    thM = 0.5 * (th[i] + th[j])
    nDotT = cEff[i] * tx + sEff[i] * ty
    rM = 0.5 * (rL[i] + rL[j])
    A = np.column_stack((thM * nDotT, tx, ty))
    b = rM * dpsi / dtP
    return A, b, psiM, thM

# Pairs with a constant on/off state; for pre-v2 logs, drop pairs whose
# c/s drive was still settling (first sample of each on/off window)
allI = np.arange(len(t) - 1)
same = (of[allI] == of[allI + 1]) & (np.diff(t) > 0)
settled = np.ones(len(t), dtype = bool)
if not fwV2:
    settled[1:] = of[1:] == of[:-1]
pairOK = same & settled[allI]
# Skip the pre-measurement stretch where c = s = 0 while of is True
pairOK &= ~(of[allI] & (c[allI] == 0.0) & (s[allI] == 0.0))
pOn = allI[pairOK & of[allI]]
pOff = allI[pairOK & ~of[allI]]

AOn1, bOn1, psiOn1, thOn1 = tangentialRows(pOn)
AOff1, bOff1, psiOff1, thOff1 = tangentialRows(pOff)
A1 = np.vstack((AOn1, AOff1))
b1v = np.concatenate((bOn1, bOff1))
co1 = np.linalg.lstsq(A1, b1v, rcond = None)[0]
g1, F0x, F0y = co1
print(f"Pairs: force-on {len(pOn)}, force-off {len(pOff)}")
print(f"Joint fit: g = {g1:.3g} /s (schematic predicts {G_PRED:.0f} /s -> "
      f"tau_eff = {5.0 / (10.0 * g1) * 1e3:.2f} ms)")
print(f"  drift F0 = ({F0x:+.4f}, {F0y:+.4f}) /s  -> |F0| = {np.hypot(F0x, F0y):.4f} /s, "
      f"attractor at {np.degrees(np.arctan2(F0y, F0x)):.1f} deg")
print(f"  equivalent velocity-node offset (tau = 10 ms): "
      f"({F0x * 5 * TAU_PRED * 1e3:+.2f}, {F0y * 5 * TAU_PRED * 1e3:+.2f}) mV")
print(f"  R^2 force-on  = {rSquared(bOn1, AOn1 @ co1):.3f}")
print(f"  R^2 force-off = {rSquared(bOff1, AOff1 @ co1):.3f}")

# Force-on only / force-off only fits for comparison
coOn1 = np.linalg.lstsq(AOn1, bOn1, rcond = None)[0]
print(f"Force-on only:  g = {coOn1[0]:.3g}, F0 = ({coOn1[1]:+.4f}, {coOn1[2]:+.4f}), "
      f"R^2 = {rSquared(bOn1, AOn1 @ coOn1):.3f}")
coOff1 = np.linalg.lstsq(AOff1, bOff1, rcond = None)[0]
print(f"Force-off only: g*C_OFF term g = {coOff1[0]:.3g}, F0 = ({coOff1[1]:+.4f}, "
      f"{coOff1[2]:+.4f}), R^2 = {rSquared(bOff1, AOff1 @ coOff1):.3f}")

# How big is the drift compared with the measurement force? Only count steps
# still collapsing (state > ~17 deg from both poles), since the tangential
# measurement force vanishes once the state reaches a pole.
measTerm = np.abs(g1 * AOn1[:, 0])
driftTerm = np.abs(AOn1[:, 1:] @ co1[1:])
nMag = np.hypot(cEff[pOn], sEff[pOn])
collapsing = np.abs(AOn1[:, 0]) > 0.3 * nMag * np.abs(thOn1)
print(f"Force-on steps still collapsing: {int(np.sum(collapsing))}")
print(f"  |measurement term| median {np.median(measTerm[collapsing]):.3f} /s, "
      f"|drift term| median {np.median(driftTerm[collapsing]):.3f} /s; "
      f"drift > measurement in {100 * np.mean(driftTerm[collapsing] > measTerm[collapsing]):.1f}% of them")
print(f"  (median |theta| where drift wins: "
      f"{np.median(np.abs(thOn1[collapsing & (driftTerm > measTerm)])):.3f}, "
      f"elsewhere: {np.median(np.abs(thOn1[collapsing & (driftTerm <= measTerm)])):.3f})")

# First-order trajectory prediction over each force-on window, driven by the
# measured theta(t) and r(t) (replaces the second-order integration above,
# whose fitted k is negative and unphysical)
def integrateFirstOrder(i0, i1, gP, fx, fy, nSub = 10):
    psiP = np.empty(i1 - i0 + 1)
    psiP[0] = psiL[i0]
    for k in range(i1 - i0):
        i = i0 + k
        h = (t[i + 1] - t[i]) / nSub
        p = psiP[k]
        for m in range(nSub):
            w = (m + 0.5) / nSub
            thW = (1 - w) * th[i] + w * th[i + 1]
            rW = (1 - w) * rL[i] + w * rL[i + 1]
            tx, ty = -np.sin(p), np.cos(p)
            p += h * (gP * thW * (cEff[i] * tx + sEff[i] * ty) + fx * tx + fy * ty) / rW
        psiP[k + 1] = p
    return psiP

traj1Err = []
traj1Examples = []
for ci, (start, end) in enumerate(mIndices):
    onIdx = start + np.where(of[start:end + 1])[0]
    if len(onIdx) < 4:
        continue
    i0 = onIdx[0] if fwV2 else onIdx[1]
    i1 = onIdx[-1]
    psiP = integrateFirstOrder(i0, i1, g1, F0x, F0y)
    traj1Err.append(np.degrees(wrapAngle(psiP[-1] - psiL[i1])))
    if len(traj1Examples) < 4 and ci >= 3:
        traj1Examples.append((ci, t[i0:i1 + 1] - t[i0], psiL[i0:i1 + 1], psiP))
traj1Err = np.array(traj1Err)
print(f"First-order trajectory prediction ({len(traj1Err)} force-on windows):")
print(f"  final angle error: median |.| = {np.median(np.abs(traj1Err)):.1f} deg, "
      f"within +/- 20 deg: {100 * np.mean(np.abs(traj1Err) <= 20):.1f}%")

t1Fig, t1Axes = plt.subplots(nrows = 1, ncols = 2, figsize = (9, 3.8))
t1Axes[0].hist(traj1Err, bins = 36, range = (-180, 180), color = 'tab:blue', alpha = 0.8)
t1Axes[0].axvline(0.0, color = 'k', ls = '--', lw = 1)
t1Axes[0].set_xlabel(r'Final angle error $\psi_{\mathrm{pred}}-\psi_{\mathrm{meas}}$ (deg)')
t1Axes[0].set_ylabel('Counts')
for k, (ci, tt, pm, pp) in enumerate(traj1Examples):
    col = f'C{k}'
    t1Axes[1].plot(tt, np.degrees(np.unwrap(pm)), color = col, lw = 2,
                   label = f'trial {ci}')
    t1Axes[1].plot(tt, np.degrees(np.unwrap(pp)), color = col, ls = '--', lw = 1.5)
t1Axes[1].set_xlabel('Time since force on (s)')
t1Axes[1].set_ylabel(r'$\psi$ (deg); solid = measured, dashed = model')
t1Axes[1].legend(fontsize = 8)
t1Fig.suptitle('First-order model: force-on trajectories')
t1Fig.tight_layout()
t1Fig.savefig('odeTrajFirstOrder.pdf')
t1Fig.show()

# Drift plot: force-off tangential speed vs angle, binned, with model
driftFig, driftAx = plt.subplots()
edges = np.linspace(-np.pi, np.pi, 25)
centers = 0.5 * (edges[1:] + edges[:-1])
binMean = np.full(len(centers), np.nan)
binSem = np.full(len(centers), np.nan)
for k in range(len(centers)):
    sel = (psiOff1 >= edges[k]) & (psiOff1 < edges[k + 1])
    if np.sum(sel) > 2:
        binMean[k] = np.mean(bOff1[sel])
        binSem[k] = np.std(bOff1[sel]) / np.sqrt(np.sum(sel))
driftAx.errorbar(np.degrees(centers), binMean, yerr = binSem, fmt = 'o',
                 ms = 4, capsize = 3, label = 'force-off data (binned)')
psiGrid = np.linspace(-np.pi, np.pi, 400)
thOffMean = float(np.mean(thOff1))
modelOff = (g1 * thOffMean * C_OFF * (np.cos(psiGrid) - np.sin(psiGrid))
            - F0x * np.sin(psiGrid) + F0y * np.cos(psiGrid))
driftAx.plot(np.degrees(psiGrid), modelOff, 'k-', lw = 1,
             label = r'joint fit (mean $\theta$)')
driftAx.axhline(0.0, color = 'gray', lw = 0.8)
driftAx.set_xlabel(r'State angle $\psi$ (deg)')
driftAx.set_ylabel(r'$r\,\dot\psi$ (fraction of 5 V per s)')
driftAx.set_title('Drift with measurement off')
driftAx.legend()
driftFig.tight_layout()
driftFig.savefig('driftOff.pdf')
driftFig.show()


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
