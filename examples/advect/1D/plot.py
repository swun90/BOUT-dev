#!/usr/bin/env python

from boutdata import collect
import matplotlib.pyplot as plt
from numpy import linspace,pi

n = collect("n")

nz = n.shape[3]

z = linspace(0,2*pi, nz, endpoint=False)

plt.plot(z, n[0,2,0,:])
plt.plot(z, n[-1,2,0,:], '-o')
plt.xlabel('Z')
plt.savefig("advect1d.pdf")
plt.show()

