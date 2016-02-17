#Time-stamp: <ssafonova.2.2.py on Tuesday, 16 February, 2016 at 18:08:33 MST (P305)>
#Sasha Safonova

#Approximates the value of the Stefan-Boltzmann constant using Romberg integration.

import numpy as np
#The script liked to throw a warning about np.exp at high values of x; this snip pet suppresses it.
import warnings; warnings.filterwarnings('ignore')

#define the romberg method.
def rom (a, b, f, error):
    romb = [[0.5 * (b-a) * (f(a) + f(b))]]
    n=1
    while True:
        h = float(b-a) / (2**n)
        romb.append((n+1)*[None])
        everything=0
        for element in range(1, 2**(n-1)+1):
            everything +=f(a+(2*element-1)*h)
        romb[n][0] = 0.5 * romb[n-1][0] + h*everything
        for i in range(1, n+1):
            romb[n][i] = romb[n][i-1] + (romb[n][i-1] - romb[n-1][i-1])/(4**i - 1)
        if np.abs(romb[n][n-1] - romb[n][n]) < error:
            return romb[n][n]
        n+=1

#define the integrand that romberg should use.
def integrand(x):
    B = ((x**3)/(np.exp(x)-1))
    return B

#declare the well-established constants
k = 1.38064852e-23
planck=6.62e-34

#Peform the integration from 0.1 to 10^5 to avoid singularities on either end of the range.
sigma = (np.pi*2*(k**4)/((planck**3)*(3e8**2)))*rom(1e-1, 1e5, integrand, 1e-11)
print ("The value of the Stefan-Bolzmann constant is: "+str(sigma))
