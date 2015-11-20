#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.optimize import minimize
from utils import fun_moog, fun_moog_fortran
import matplotlib.pyplot as plt



Teff_rng = np.arange(4800, 5000, 50)
logg_rng = np.arange(4.0, 4.3, 0.1)


## Python version
Teffs, loggs = np.meshgrid(Teff_rng, logg_rng)
z = np.zeros(Teffs.T.shape)
for i, Teff in enumerate(Teffs[0, :]):
    for j, logg in enumerate(loggs[:, 0]):
        print 'Python --\tTeff: %i/%i\tlogg: %.2f/%.2f' % (Teff, Teff_rng.max(), logg, logg_rng.max())
        z[i, j] = fun_moog((Teff, logg, 0.0, 1.0))[0]

z = z.T
z1 = z.copy()
z /= z.max()

## FORTRAN version
# Teffs, loggs = np.meshgrid(Teff_rng, logg_rng)
#
# zf = np.zeros(loggs.T.shape)
# for i, Teff in enumerate(Teffs[0, :]):
#     for j, logg in enumerate(loggs[:, 0]):
#         print 'Fortran --\tTeff: %i/%i\tlogg: %.2f/%.2f' % (Teff, Teff_rng.max(), logg, logg_rng.max())
#         zf[i, j] = fun_moog_fortran((Teff, logg, 0.0, 1.0))[0]
#
# zf = zf.T
# z2 = zf.copy()
# zf /= zf.max()

# res = abs(z1-z2)

plt.figure()
plt.contourf(Teffs, loggs, z, levels=np.linspace(0, 1, 100), cmap=plt.cm.spectral)
plt.colorbar()
plt.title('Python')
# plt.plot(5834, 4.28, 'or')


# plt.figure()
# plt.contourf(Teffs, loggs, zf, levels=np.linspace(0, 1, 100), cmap=plt.cm.spectral)
# plt.colorbar()
# plt.title('Fortran')
# # # plt.plot(5834, 4.28, 'or')
# #
# plt.figure()
# plt.contourf(Teffs, loggs, res, levels=np.linspace(0, 1, 100), cmap=plt.cm.spectral)
# plt.colorbar()
# plt.title('Difference')
# # plt.plot(5834, 4.28, 'or')

plt.show()
