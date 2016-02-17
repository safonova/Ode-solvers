#Time-stamp: <ps2_1.py on Tuesday, 9 February, 2016 at 12:11:35 MST (P305)>
#Sasha Safonova

#Integrates \int_{0}^{1}x*e^{-x}dx using tiny trapezoids
#Area of a trapezoid:
#A=h*(a+b)/2

import numpy as np

def integrand(x):
    return x*np.exp(-x)

def traprule(f, low, high, n):
    h = (high - low) / n
    I = 0.5 * (f(low) - f(high)) * h
    for i in range(1, n-2):
        I += h * f(low + i*h)
    return I

answer = traprule(integrand, 0, 1, 1000240)
print ("The integral of x*e^-x from 0 to 1 is approximately "+str(answer))
