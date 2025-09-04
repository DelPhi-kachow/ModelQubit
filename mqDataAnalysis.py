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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

#Setup variables
fname = 'modelQubit_09-03-2025_182741.log'

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
        print(row)
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
count = 0
mIndices = []
incycle = False
startdex = 0
enddex = 0
#Iterate through on/off
for i, elem in enumerate(of):
    #print(i, elem)
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
            count += 1
            mIndices.append( (startdex, enddex) )
            
            #Begin the next cycle
            startdex = i

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

#Try everything out for a specific measurement
#startdex = mIndices[81][0]
#enddex = mIndices[81][1]

#for i in range(startdex, enddex):
    #print(of[i], c[i], s[i], x[i], y[i], updown(i))
    
mmts = np.array([updown(i[1]) for i in mIndices])
#print(mmts)

#I'm seeing a couple 0s, which is not good
#Let's see if we can figure out why
#I need to locate and plot one of these
'''
After doing the plotting, it seems clear that (x,y) was on its way to the 
expected location, but the measurement ended too early. To see how often this
happened, I can compute the expected measurement from theta and see how often
the test had enough time.
'''

#Time average theta during a trial
theets = np.array([
    np.average(th[mIndices[i][0] : mIndices[i][1]]) 
    for i in range(len(mIndices))
    ])
#The measurement is expected to converge to the sign of theta
expectation = np.sign(theets)
#print(expectation)

#How often was that the case?
deviation = mmts - expectation
#print(deviation)
failFreq = np.sum(deviation / 2.) / len(deviation)
#print(failFreq)

#This is over a 30% failure rate!
#I need to double check some of the identified failures in animation
#It does indeed look like many tests did not have enough time under measurement
#I'll need to redo the tests with more measurement time
#Okay, I added more time, but now the failure rate is still about 13%
#Let's double check that the current failures are still the same sort of thing
fails = np.nonzero(deviation)
#print(fails)
#Animating a few of these, and it looks like now we're seeing the measurement
# go on long enough that the sign of theta is changing. I think this is fine now.

#Bin the relative angles and compute measured probabilities
binCount = 10
binStep = (np.pi / 2.0) / (binCount - 1)
inBins = np.zeros(binCount) #Storing how many measurements in bin
bins = np.zeros(binCount) #Storing net measurement of bin
for i in range(len(mIndices)):
    angle = rangles[i]
    #Which bin does the angle go in?
    binDex = int(np.round( angle / binStep ))
    #Count trial
    inBins[binDex] += 1
    #Count result
    bins[binDex] += (mmts[i] + 1.) / 2. 

#Compute probabilities of measuring up
probs = np.array([bins[i] / inBins[i] for i in range(binCount)])
print(inBins)
print(probs)


'''
PLOT MEASURED VS IDEAL QUBIT
'''

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

#Plot an xy point
xy, = ax.plot([x[1000]], [y[1000]], marker = 'o')

#Draw the measurement axis
axup, = ax.plot([0., c[1000]], [0., s[1000]], color = 'orange')
axdown, = ax.plot([0., -c[1000]], [0., -s[1000]], color = 'orange', ls = '--')
axtip, = ax.plot(c[1000], s[1000], color = 'orange', marker = 'o', mfc = 'white')

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
'''

'''
#Animation of full trial
for i in range(len(mIndices)):
    print('Making video for trial ' + str(i))
    
    start = mIndices[i][0]
    end = mIndices[i][1]
    video = ani.FuncAnimation(fig, step,
                              repeat = True, 
                              frames = range(start - 40, end + 40),
                              interval = 50,
                              blit = True)
    #fig.show()
    
    #To save the animation using Pillow as a gif
    writer = ani.PillowWriter(fps=20)
    ffname = 'test_' + str(i) + '.gif'
    video.save(ffname, writer=writer)
'''

'''
#Animation of a single test, where neither up or down was measured
start = mIndices[740][0]
end = mIndices[740][1]
video = ani.FuncAnimation(fig, step,
                          repeat = True, 
                          frames = range(start - 40, end + 40),
                          interval = 50,
                          blit = True)
#fig.show()

#To save the animation using Pillow as a gif
writer = ani.PillowWriter(fps=20)
video.save('weirdTest5.gif', writer=writer)

#It's clear that it was on its way to the expected state, but the measurement
# period ended too early.
'''