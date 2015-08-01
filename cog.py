#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    d = np.loadtxt(sys.argv[1], skiprows=1)
else:
    d = np.loadtxt('linelist.moog', skiprows=1)

d = d[d[:,1] == 26.0]

rw = np.log(d[:,-1]/d[:,0])
def A(teff):
    return d[:,3] - 5040/teff * d[:,2]

teffs = np.arange(3000, 10000, 10)
std = np.zeros(teffs.shape)

for i, teff in enumerate(teffs):
    ai = A(teff)
    p = np.polyfit(ai, rw, 2)
    res = rw - np.poly1d(p)(ai)
    std[i] = np.std(res)
    # if not (i % 5) and std[i] < 0.1:
    #     plt.subplot(211)
    #     plt.plot(ai, rw, 'ok')
    #     plt.plot(sorted(ai), np.poly1d(p)(sorted(ai)), '-r')
    #     plt.subplot(212)
    #     plt.plot(ai, res, '.k')
    #     plt.show()

T = teffs[np.argmin(std)]
print 'Best Teff from CoG: %iK' % T
plt.plot(teffs, std, '-k', lw=2)
plt.title('Curve of growth')
plt.xlabel('Teff [K]')
plt.ylabel('Residuals: CoG - linear fit')

plt.figure()
plt.plot(A(5777), rw, 'ok')
plt.title('5777')

plt.figure()
plt.plot(A(T), rw, 'ok')
plt.title('%i' % T)
plt.show()