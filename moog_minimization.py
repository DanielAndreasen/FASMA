#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from runmoog import fun_moog_fortran, fun_moog


def print_format(x):
    print '%i, %.2f, %.2f, %.2f' % (x[0], x[1], x[2], x[3])


def minimize(x0, func, bounds=None, fix_logg=False, eps=0.0001, iteration=100):
    """
    Some doc
    """

    # Initial step
    # Teff, logg, [Fe/H], vt
    step = (500, 0.10, 0.10, 0.10)

    res, slopeEP, slopeRW, abundances = func(x0)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)

    N = 0
    while (res > eps) and (N < iteration):
        while (abs(slopeEP) > 0.001):
            # For Teff
            print 'TEMP'
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**2
            parameters[0] += step_i
            # parameters[2] = abundances[0] - 7.47
            print_format(parameters)
            print N, s, slopeEP, res, '\n'
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

        while (abs(Abdiff) > 0.001) and not fix_logg:
            # For logg
            print 'LOGG'
            s = -np.sign(Abdiff)
            step_i = s * step[1]/abs(np.log(abs(Abdiff)+0.0005))**2
            parameters[1] += step_i
            # parameters[2] = abundances[0] - 7.47
            print_format(parameters)
            print N, s, Abdiff, res, '\n'
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

        while (abs(slopeRW) > 0.001):
            # For micro turbulence
            print 'VT'
            s = np.sign(slopeRW)
            step_i = s * step[3]/abs(np.log(abs(slopeRW)+0.0005))**2
            parameters[3] += step_i
            # parameters[2] = abundances[0] - 7.47
            print_format(parameters)
            print N, s, slopeRW, res, '\n'
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

    return parameters


x0 = (5700, 4.44, 0.0, 1.0)
minimize(x0, fun_moog_fortran, fix_logg=True)
