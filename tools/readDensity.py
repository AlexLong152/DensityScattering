# -*- coding: utf-8 -*-

"""
@author: alexl
"""
import numpy as np
import h5py
filename = r'/home/alex/OneDrive/densities-3He/2Ndensities/twoBody-3He-theta=59.98-omega=131.857-chsfr-nnlo-lambda=2.0-hash=4dc5aaae23f81daea252bd01a7ad1b53de75e0e173eef9ecdc651b8284a4d8ee-row=106.h5'


# f = h5py.File(filename, "r")
with h5py.File(filename, "r") as f:
    print("Keys: %s" % f.keys())

    keys = list(f.keys())
    print("keys[0]=",keys[0])
    rhoGroup = f[keys[0]]
    rhoKeys = list(rhoGroup.keys())
    p12p = f['p12p'][()]
    rho = rhoGroup['RHO'][()]
print(np.shape(rho))
for i in range(len(rho)):
    print(rho[0,i,i])
# [()] dumps this into numpy array that can be used later

# f.close()
