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

#For each measurement, what was the relative angle between (x, y) and (c, s)?

#Did we measure up or down on this axis?

#Bin the relative angles and compute measured probabilities

'''
PLOT MEASURED VS IDEAL QUBIT
'''

'''
ANIMATIONS
'''

