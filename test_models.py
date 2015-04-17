#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import seaborn as sns
sns.set_style('dark')
sns.set_context('talk')
import matplotlib.pyplot as plt
from model_interpolation import _unpack_model, read_model
import sys

py = sys.argv[1]
fo = sys.argv[2]


d1 = np.loadtxt(py)
d2 = np.genfromtxt(fo, skip_header=3, skip_footer=7)
# d2 = np.genfromtxt(fo, skip_header=3, skip_footer=7)


#plt.plot(d1[:, 0], 100*(d1[:, 1]-d2[:,1])/d1[:,1], 'o')
#plt.hlines(np.median(100*(d1[:, 1]-d2[:,1])/d1[:,1]), d1[0, 0], d1[-1, 0], 'r')
plt.plot(d2[:, 0], d2[:, 1], 'o',d1[:, 0], d1[:, 1], 'o')
plt.title("Temperature")

plt.figure()

#plt.plot(d1[:, 0], 100*(d1[:, 2]-d2[:,2])/d1[:,2], 'o') 
#plt.hlines(np.median(100*(d1[:, 2]-d2[:,2])/d1[:,2]), d1[0, 0], d1[-1, 0], 'r')
plt.plot(d2[:, 0], d2[:, 2], 'o',d1[:, 0], d1[:, 2], 'o')
plt.title("Presssure")


plt.figure()
#plt.plot(d1[:, 0], 100*(d1[:, 3]-d2[:,3])/d1[:,3], 'o')
#plt.hlines(np.median(100*(d1[:, 3]-d2[:,3])/d1[:,3]), d1[0, 0], d1[-1, 0], 'r')
plt.plot(d2[:, 0], d2[:, 3], 'o',d1[:, 0], d1[:, 3], 'o')
plt.title("Numb of Electrons")

plt.figure()
#plt.plot(d1[:, 0], 100*(d1[:, 4]-d2[:,4])/d1[:,4], 'o')
#plt.hlines(np.median(100*(d1[:, 4]-d2[:,4])/d1[:,4]), d1[0, 0], d1[-1, 0], 'r')
plt.plot(d2[:, 0], d2[:, 4], 'o',d1[:, 0], d1[:, 4], 'o')
plt.title("Rossland absorption")

plt.figure()
#plt.plot(d1[:, 0], 100*(d1[:, 5]-d2[:,5])/d1[:,5], 'o')
#plt.hlines(np.median(100*(d1[:, 5]-d2[:,5])/d1[:,5]), d1[0, 0], d1[-1, 0], 'r')

plt.plot(d2[:, 0], d2[:, 5], 'o',d1[:, 0], d1[:, 5], 'o')
plt.title("Accrad?")

plt.show()




#plt.plot(d2[:,0],d2[:,0]-d1[:,0], 'o')
#plt.show()
