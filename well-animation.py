#!/usr/bin/python3

# Program that calculates the evolution of a particle initially
# confined to the right half of an infinite square well.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

from math import *

hbar = 1.05e-34 # J*s
a = 1e-9 # m
m = 9.11e-31 # kg
nmax = 1000
dt = 1e-19 # s
nframes = 1000
npoints = 100

def E(n):
    return hbar**2*pi**2*n**2/(2*m*a**2)

def c(n):
    return sqrt(2/a)*(-a/(n*pi))*(cos(n*pi)-cos(n*pi/2))

def psin(n,x):
    return sqrt(2/a)*np.sin(n*pi*x/a)

def psi(x,t):
    total=0
    for n in range(1,nmax+1):
        total=total+c(n)*psin(n,x)*np.exp(-1j*E(n)*t/hbar)
    return total

def psisquared(x,t):
    return abs(psi(x,t))**2


fig, ax = plt.subplots()

x = np.linspace(0, a, npoints)

t = 0

line2 = ax.plot(x, psisquared(x,t))[0]
ax.set(xlabel='x [m]', ylabel='Prob density')
# plt.show() # use to test

def update(frame):
    # print('Frame ',frame)
    t = dt*frame
    line2.set_ydata(psisquared(x,t))
    return line2,


ani = animation.FuncAnimation(fig=fig, func=update, frames=nframes, interval=.1)
writergif = animation.PillowWriter(fps=30)
ani.save('well-animation.gif', writer = writergif)
plt.show()
