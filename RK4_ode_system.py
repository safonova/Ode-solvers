#Time-stamp: <ssafonova.2.6.py on Tuesday, 16 February, 2016 at 23:27:36 MST (P305)>
'''Sasha Safonova
PHYS305

Use Runge Kutta-4 to solve a system of two first-order ode's stemming from a second-order equation:
y''(t)+\omega^2*y(t)=0,
where \omega=1

initial conditions:
y(0)=0
y'(0)=1
y''(0)=0'''

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm

#rhs particular to our case.
def rhs(x, y, z):
    x = x
    f = z
    g = -y
    return x, f, g

#finds RK4 solutions for a system of 2 equations
def RK4 (init, step, derivs):
    xi, yi, zi = init
    x0, k0, l0 = derivs(xi, yi, zi)
    k0 = k0*step
    l0 = l0*step
    x1, k1, l1 = derivs(xi + 0.5*step, yi+0.5*k0, zi+0.5*l0)
    k1 = k1*step
    l1 = l1*step
    x2, k2, l2 = derivs(xi + 0.5*step, yi+0.5*k1, zi+0.5*l1)
    k2 = k2*step
    l2 = l2*step
    x3, k3, l3 = derivs(xi + step, yi+k2, zi+l2)
    k3 = k3*step
    l3 = l3*step
    ynext = yi + (k0 + 2*k1 + 2*k2 + k3)/6
    znext = zi + (l0 + 2*l1 + 2*l2 + l3)/6
    return [xi, ynext, znext]

#method for finding RK4 solutions over the desired range of time
def execute(seconds, step):
    time=0
    init=[ 0, 0, 1]
    a=np.linspace(time, seconds, ((seconds-time)/step))
    results=[init]
    for n in a:
        init[0]=n
        init = RK4(init, step, rhs)
        results.append(init)
    return np.array(results)

#the pen-and-paper solution
def analytic(x):
    return np.sin(x)

#find one set of solutions and plot RK4 and the analytic solutions
large = execute(16*np.pi, 0.1)

plt.figure()
z=range(len(large[:,0]))
a=np.linspace(0,(16*np.pi))
plt.scatter(large[:,0], large[:,1], c=z, cmap=cm.get_cmap('rainbow'), label="RK4", edgecolor='none')
plt.plot(a, analytic(a), color='skyblue', linewidth=5, alpha=0.4, label="Analytic")
plt.title("$y''(t)+\omega^2*y(t)=0$")
plt.legend(loc='best')
plt.show()
