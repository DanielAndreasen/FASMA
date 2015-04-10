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


d1 = np.genfromtxt(py, skip_header=3, skip_footer=7)
d2 = np.genfromtxt(fo, skip_header=3, skip_footer=7)

plt.plot(d1[:, 0], d1[:, 1], 'o')
plt.plot(d2[:, 0], d2[:, 1], 'o')
plt.show()
