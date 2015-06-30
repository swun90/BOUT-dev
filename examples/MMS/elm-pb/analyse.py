#!/usr/bin/env python

# Python script to analyse MMS test output

from boutdata import collect

from numpy import sqrt, max, abs, mean, array

import pickle

nxlist = [8, 16, 32, 64]
varlist = ["P", "Psi", "U"]

tind = 1

error_2 = {}
error_inf = {}

for nx in nxlist:
    directory = "grid%d" % nx
    
    for var in varlist:
        # Collect data
        E = collect("E_"+var, tind=[tind,tind], info=False, path=directory)
        E = E[:,2:-2, :,:]
        
        # Average error over domain
        l2 = sqrt(mean(E**2))
        linf = max(abs( E ))
    
        error_2[var] = l2
        error_inf[var] = linf

        print("%s : l-2 %f l-inf %f" % (var, l2, linf))
    
    # Save data
    with open(directory+".pkl", "wb") as output:
        pickle.dump(nx, output)
        pickle.dump(error_2, output)
        pickle.dump(error_inf, output)
    
