#!/usr/bin/env python

from boutdata import collect

import matplotlib.pyplot as plt

n = collect("n")
#psi = collect("psi")

cmap=None
#cmap = cm.gray
plt.imshow(n[-1,:,0,:], interpolation='bilinear', cmap=cmap, origin='lower')
#plt.contourf(n[-1,:,0,:])
plt.contour(n[0,:,0,:],levels=[0.5])

#plt.contour(psi[:,0,:])

#plt.colorbar(shrink=0.8, extend='both')

plt.xlabel("X")
plt.ylabel("Z")

plt.savefig("advect2d.pdf")

from numpy import amin, amax

print "MIN: ", amin(n[-1,:,0,:])
print "MAX: ", amax(n[-1,:,0,:])

plt.show()
