#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from runmoog import fun_moog_fortran, fun_moog


def print_format(x):
    print '%i, %.2f, %.2f, %.2f' % (x[0], x[1], x[2], x[3])


def minimize(x0, func, bounds=None,
             fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
             eps=1e-7, iteration=100):
    """
    Some doc
    """

    # TODO: Implement the bound. e.g. if we hit upper Teff, raise a flag and
    # set the Teff=upper bound

    # Step size
    # Teff, logg, vt
    step = (500, 0.50, 0.50)

    res, slopeEP, slopeRW, abundances = func(x0)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)

    N = 0
    while (res > eps) and (N < iteration):
        while (abs(slopeEP) > 0.001) and not fix_teff:
            # For Teff
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**3
            step_i = s*1 if abs(step_i) < 1 else step_i
            parameters[0] += step_i
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

        while (abs(Abdiff) > 0.001) and not fix_logg:
            # For logg
            s = -np.sign(Abdiff)
            step_i = s * step[1]/abs(np.log(abs(Abdiff)+0.0005))**3
            step_i = s*0.01 if abs(step_i) < 0.01 else step_i
            parameters[1] += step_i
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

        #       Input metal...   FeI abund.
        while (parameters[2] != abundances[0]-7.47) and not fix_feh:
            # For metalicity
            parameters[2] = abundances[0]-7.47
            res, slopeEP, slopeRW, abundances = func(parameters)
            print_format(parameters)
        N += 1

        while (abs(slopeRW) > 0.001) and not fix_vt:
            # For micro turbulence
            s = np.sign(slopeRW)
            step_i = s * step[2]/abs(np.log(abs(slopeRW)+0.0005))**3
            step_i = s*0.01 if abs(step_i) < 0.01 else step_i
            parameters[3] += step_i
            parameters[3] = 0.00 if parameters[3] < 0 else parameters[3]
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
        N += 1

    return parameters


x0 = (5777, 4.44, 0.0, 1.00)
# minimize(x0, fun_moog_fortran)
minimize(x0, fun_moog)
