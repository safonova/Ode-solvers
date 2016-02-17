#Time-stamp: <ssafonova.2.4.py on Monday, 15 February, 2016 at 21:38:31 MST (P305)>
#Sasha Safonova

#integrate df/dt=-f/7.7, f(0)=1 using the forward Euler method
#Use the formula f(t+h)=f(t)+c*h*f(t)

import numpy as np
from matplotlib import pyplot as plt

#Declare the constant coefficient of our current ode's rhs.
def rhs ():
    c= (-1/7.7)
    return c

#method for implementing the Euler from line 5
def forwardeuler (x, y, rhs, step):
    f_t_h = y + step*rhs()*y
    return f_t_h

#declare the "correct" analytic solution
def analytic(x):
    y = np.exp(-x/7.7)
    return y

#
def data(days, step):
    x=0
    y=1
    a=np.linspace(x, days, (days/step))
    results = []
    for n in a:
        results.append([x,y])
        result = forwardeuler(x, y, rhs, step)
        y=result
        x+=step
    results=np.array(results)
    print ("f("+str(x)+") = "+str(y))
    return results

#gather data for different step sizes
large = data(5, 1)
medium = data(5, 0.1)
small = data(5, 1e-2)

#plot the datasets in different alpha values to make the different apparent
plt.figure()
plt.scatter(large[:,0], large[:,1], c='b',alpha=0.6, label="1 step")
plt.scatter(medium[:,0], medium[:,1], c='g',alpha=0.4, label="0.1 step")
plt.scatter(small[:,0], small[:,1], c='r',alpha=0.2,edgecolor='none', label="1e-2 step")
plt.plot(small[:1],analytic(small[:1]), c='c', label="analytic")
plt.title("Decay function")
plt.legend(loc="best")
plt.show()
