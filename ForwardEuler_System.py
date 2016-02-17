#Time-stamp: <ssafonova.2.5.py on Tuesday, 16 February, 2016 at 23:30:54 MST (P305)>
'''Sasha Safonova
PHYS 305
Solves a system of ode using the forward Euler method
Uses the formula f(t+h)=f(t)+c*h*f(t)'''

import numpy as np
from matplotlib import pyplot as plt

#Declare a system of odes.
def system (vector):
    a = (-1/7.7)*vector[0]
    b = (1/7.7)*vector[0] -(1/111)*vector[1]
    return a, b

#method for implementing the Euler from line 5 using a vector for the ode system
def forwardeuler (yvector, system, step):
    newvector=[]
    euler=system(yvector)
    for k in range(len(yvector)):
        l=yvector[k]+step*euler[k]
        newvector.append(l)
    return newvector

#generates step-by-step solution array
def data(days, step):
    #x=[0]*3
    x=0
    y=[1,0]
    a=np.linspace(x, days, ((days-x)/step))
    results = []
    for n in a:
        ni, co= forwardeuler(y, system, step)
        fe=1-ni-co
        y=[ni, co]
        results.append([n,ni,co,fe])
    results=np.array(results)
    return results

#gather data for different step sizes
large = data(5, 1)
medium = data(5, 0.1)
small = data(5, 1e-2)

#plot the datasets in different alpha values to make the different apparent
fig=plt.figure()
axni = fig.add_subplot(221)
axni.set_title("Ni vs time (days)")
axco = fig.add_subplot(222)
axco.set_title("Co vs time (days)")
axfe = fig.add_subplot(212)
axfe.set_title("Fe vs time (days)")

axni.scatter(large[:,0], large[:,1], c='skyblue',alpha=0.6, label="Ni, 1 step")
axco.scatter(large[:,0], large[:,2], c='peru',alpha=0.6, label="Co, 1 step")
axfe.scatter(large[:,0], large[:,3], c='lightslategray',alpha=0.6, label="Fe, 1 step")
axni.scatter(medium[:,0], medium[:,1], c='royalblue',alpha=0.4, label="0.1 step")
axco.scatter(medium[:,0], medium[:,2], c='darkgoldenrod',alpha=0.4, label="0.1 step")
axfe.scatter(medium[:,0], medium[:,3], c='lightgrey',alpha=0.4, label="0.1 step")
axni.scatter(small[:,0], small[:,1], c='cornflowerblue',alpha=0.1, label="1e-2 step", edgecolor='none')
axco.scatter(small[:,0], small[:,2], c='brown',alpha=0.1, label="1e-2 step", edgecolor='none')
axfe.scatter(small[:,0], small[:,3], c='silver',alpha=0.1, label="1e-2 step", edgecolor='none')


axni.legend(loc=3, prop={'size':8})
axco.legend(loc=4, prop={'size':8})
axfe.legend(loc="lower center", prop={'size':8})
plt.show()
