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


def plot_absolute(d1, d2, col, title=False):
    plt.figure()
    plt.plot(d1[:, 0], d1[:, col], '-o', label='Python')
    plt.plot(d2[:, 0], d2[:, col], '-o', label='FORTRAN')
    if title:
        plt.title(title)
    plt.legend()


def plot_relative(d1, d2, col, title=False):
    plt.figure()
    plt.plot(d1[:, 0], 100*(d1[:, col]-d2[:, col])/d1[:, col], '-o')
    plt.hlines(np.median(100*(d1[:, col]-d2[:, col])/d1[:, col]), d1[0, 0], d1[-1, 0], 'r')
    if title:
        plt.title(title)


py = sys.argv[1]
fo = sys.argv[2]


d1 = np.loadtxt(py)
d2 = np.genfromtxt(fo, skip_header=3, skip_footer=7)
# d2 = np.genfromtxt(fo, skip_header=3, skip_footer=7)

titles = 'Temperature,Presssure,Number of electrons, Rossland absorption,Accrad'.split(',')

for col, title in enumerate(titles):
    # plot_absolute(d1, d2, col+1, title)
    plot_relative(d1, d2, col+1, title)
plt.show()
