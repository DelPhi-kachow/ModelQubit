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
ANIMATIONS
'''

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

#Animation of full trial
for i in range(len(mIndices)):
    ffname = 'test_' + str(i) + '.gif'
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
'''

#Optional: animate a specific inconclusive trial if it exists in this log
#It's clear that it was on its way to the expected state, but the measurement
# period ended too early.
'''
weirdDex = 740
if weirdDex < len(mIndices) and not os.path.exists('weirdTest5.gif'):
    start = mIndices[weirdDex][0]
    end = mIndices[weirdDex][1]
    video = ani.FuncAnimation(fig, step,
                              repeat = True, 
                              frames = trial_frames(start, end),
                              interval = 50,
                              blit = True)
    writer = ani.PillowWriter(fps=20)
    video.save('weirdTest5.gif', writer=writer)
'''
