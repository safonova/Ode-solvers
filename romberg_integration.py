#Time-stamp: <ps2_3.py on Sunday, 14 February, 2016 at 17:30:47 MST (P305)>
#Sasha Safonova

#Approximates the value of e^(-x^2) over all space.

import numpy as np

#define the romberg method (as written in ps2_2.py)
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
    B = np.exp(-x**2)
    return B


'''Peform the integration from 0 to 10^4 and multiply by 2 to exploit symmetry in the solution.
Then compare with the value that Lord Kelvin thought was so obvious.'''
answer = 2*rom(0, 1e4, integrand, 1e-7)
print ("The value is: "+str(answer))
print ("Pi over 2 is: "+str(np.pi/2))
