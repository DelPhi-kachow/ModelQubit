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

#Setup variables
fname = 'modelQubit_08-25-2025_230700.log'

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
            enddex = i
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

rangles = np.array([rangle(dex[0]) for dex in mIndices])

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
startdex = mIndices[1][0]
enddex = mIndices[1][1]

for i in range(startdex, enddex):
    print(of[i], c[i], s[i], x[i], y[i], updown(i))
    
mmts = np.array([updown(i[1]-1) for i in mIndices])
print(mmts)

#I'm seeing a couple 0s, which is not good
#Let's see if we can figure out why
#I need to locate and plot one of these


#Bin the relative angles and compute measured probabilities

'''
PLOT MEASURED VS IDEAL QUBIT
'''

'''
ANIMATIONS
'''

