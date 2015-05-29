#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.optimize import minimize
from runmoog import fun_moog, fun_moog_fortran
import matplotlib.pyplot as plt



# Teff_rng = np.arange(4500, 6500, 500)
# logg_rng = np.arange(4.0, 5.0, 0.1)

## Python version
# Teffs, loggs = np.meshgrid(Teff_rng, logg_rng)
# z = np.zeros(Teffs.T.shape)
# for i, Teff in enumerate(Teffs[0, :]):
#     for j, logg in enumerate(loggs[:, 0]):
#         print Teff, logg
#         z[i, j] = fun_moog((Teff/1000, logg, 0.0, 1.0))

# z = z.T
# plt.contourf(Teffs, loggs, z)
# plt.colorbar()
# plt.plot(5777, 4.44, 'ok')



## FORTRAN version
# Teff_rng = np.arange(4000, 8000, 50)
logg_rng = np.arange(3.1, 4.8, 0.1)
feh_rng = np.arange(-2.0, 1.5, 0.1)
# vt_rng = np.arange(0.1, 4.0, 0.1)
loggs, fehs = np.meshgrid(logg_rng, feh_rng)
zf = np.zeros(loggs.T.shape)
for i, logg in enumerate(loggs[0, :]):
    for j, feh in enumerate(fehs[:, 0]):
        print logg, feh
        zf[i, j] = fun_moog_fortran((5777, logg, feh, 1.0))

zf = zf.T
zf /= zf[~np.isnan(zf)].max()
plt.figure()
plt.contourf(loggs, fehs, zf, levels=np.linspace(0, 1, 100),
        cmap=plt.cm.spectral)
plt.colorbar()
plt.plot(4.44, 0.0, 'or')
plt.show()
