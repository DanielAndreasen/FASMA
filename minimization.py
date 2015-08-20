#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os


def print_format(x):
    """Print the stellar atmospheric parameters in a nice format"""
    print '%i, %.2f, %.2f, %.2f' % (x[0], x[1], x[2], x[3])


def save_iteration(parameters):
    """Save the current iteration"""
    data = ','.join(map(str, parameters))
    os.system('echo %s >> minimization_profile.dat' % data)


def check_bounds(parameter, bounds, i):
    """
    Function which checks if parameters are within bounds.
    Input - parameter: what we want to check; bounds: ze bounds;
    i: the index of the bounds we want to check"""
    if parameter < bounds[i-1]:
        parameter = bounds[i-1]
    elif parameter > bounds[i]:
        parameter = bounds[i]
    return parameter


def check_convergence(RW, EP, Abdiff, fe_input, fe,
                      fix_teff, fix_logg, fix_vt, fix_feh):
    """
    Check convergence criteria
    """
    EP = 0.00 if fix_teff else EP
    RW = 0.00 if fix_vt else RW
    Abdiff = 0.00 if fix_logg else Abdiff
    fe = fe_input+7.47 if fix_feh else fe

    cond1 = abs(RW) <= 0.001
    cond2 = abs(Abdiff) <= 0.01
    cond3 = abs(EP) <= 0.001
    cond4 = fe_input == fe-7.47
    return cond1 and cond2 and cond3 and cond4


def minimize(x0, func, bounds="kurucz95",
             fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
             iteration=25):
    """
    Sane minimization like a normal human being would do it.
    """

    # Teff, logg, vt
    step = (500, 0.50, 0.50)
    if bounds.lower() == "kurucz95":
        bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]

    res, slopeEP, slopeRW, abundances = func(x0)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)
    try:
        os.remove('minimization_profile.dat')
    except OSError:
        pass

    N = 0
    while N < iteration:
        N1 = 0
        cycle = [parameters[0]]
        while (abs(slopeEP) > 0.001) and not fix_teff and N1 < 15:
            # For Teff
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**3
            if abs(step_i) < 1:
                step_i = s*1  # Minimum 1K
            elif abs(step_i) > 17625:
                step_i = s*17625  # Maximum 17625K (half the interval?)
            parameters[0] += step_i
            parameters[0] = check_bounds(parameters[0], bounds, 1)
            if parameters[0] in cycle:
                break
            else:
                cycle.append(parameters[0])
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
            N1 += 1
            save_iteration(parameters)
        N += 1
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0], fix_teff, fix_logg, fix_vt, fix_feh):
            break

        N4 = 0
        cycle = [parameters[3]]
        while (abs(slopeRW) > 0.001) and not fix_vt and N4 < 15:
            # For micro turbulence
            s = np.sign(slopeRW)
            step_i = s * step[2]/abs(np.log(abs(slopeRW)+0.0005))**3
            if abs(step_i) < 0.01:
                step_i = s*0.01  # Minimum vt step
            elif abs(step_i) > 2.5:
                step_i = s*2.5  # Maximum vt step
            else:
                step_i = step_i
            parameters[3] += step_i
            parameters[3] = 0.00 if parameters[3] < 0 else parameters[3]
            parameters[3] = check_bounds(parameters[3], bounds, 7)
            if parameters[3] in cycle:
                break
            else:
                cycle.append(parameters[3])
            print_format(parameters)
            tmp = slopeRW
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
            N4 += 1
            save_iteration(parameters)
            if tmp == slopeRW:
                break

        N2 = 0
        cycle = [parameters[1]]
        while (abs(Abdiff) > 0.01) and not fix_logg and N2 < 15:
            # For logg
            s = -np.sign(Abdiff)
            step_i = s * step[1]/abs(np.log(abs(Abdiff)+0.0005))**3
            if abs(step_i) < 0.01:
                step_i = s*0.01  # Minimum logg step
            elif abs(step_i) > 2.5:
                step_i = s*2.5  # Maximum logg step
            parameters[1] += step_i
            # checks bounds of logg
            parameters[1] = check_bounds(parameters[1], bounds, 3)
            if parameters[1] in cycle:
                break
            else:
                cycle.append(parameters[1])
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters)
            Abdiff = np.diff(abundances)[0]
            N2 += 1
            save_iteration(parameters)
        N += 1
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0], fix_teff, fix_logg, fix_vt, fix_feh):
            break

        #       Input metal...   FeI abund.
        N3 = 0
        while (parameters[2] != abundances[0]-7.47) and not fix_feh and N3 < 5:
            # For metalicity
            parameters[2] = abundances[0]-7.47
            parameters[2] = check_bounds(parameters[2], bounds, 5)
            res, slopeEP, slopeRW, abundances = func(parameters)
            print_format(parameters)
            N3 += 1
            save_iteration(parameters)
        N += 1
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0], fix_teff, fix_logg, fix_vt, fix_feh):
            break

        N += 1
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0], fix_teff, fix_logg, fix_vt, fix_feh):
            break

    print 'Stopped at %i iterations' % N
    converged = check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0], fix_teff, fix_logg, fix_vt, fix_feh)
    parameters[0] = int(parameters[0])
    parameters[1] = round(parameters[1], 2)
    parameters[2] = round(parameters[2], 2)
    parameters[3] = round(parameters[3], 2)

    return parameters, converged
