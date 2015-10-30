#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from copy import copy


def print_format(x):
    """Print the stellar atmospheric parameters in a nice format"""
    print '%i %.2f %.2f %.2f' % (x[0], x[1], x[2], x[3])


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


def check_convergence(RW, EP, Abdiff, fe_input, fe):
    """
    Check convergence criteria
    """
    EP = 0.00 if f_teff else EP
    RW = 0.00 if f_vt else RW
    Abdiff = 0.00 if f_logg else Abdiff
    fe = fe_input+7.47 if f_feh else fe

    cond1 = abs(RW) <= RWcriteria
    cond2 = abs(Abdiff) <= ABdiffcriteria
    cond3 = abs(EP) <= EPcriteria
    cond4 = round(fe_input, 2) == round(fe-7.47, 2)
    return cond1 and cond2 and cond3 and cond4


def minimize(x0, func, bounds="kurucz95", weights='null',
             fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
             iteration=160, EPcrit=0.002, RWcrit=0.003, ABdiffcrit=0.01):
    """
    Sane minimization like a normal human being would do it.
    """

    global EPcriteria, RWcriteria, ABdiffcriteria
    global f_teff, f_logg, f_vt, f_feh
    EPcriteria, RWcriteria, ABdiffcriteria = EPcrit, RWcrit, ABdiffcrit
    f_teff, f_logg, f_feh, f_vt = fix_teff, fix_logg, fix_feh, fix_vt
    # Step size in Teff, logg, vt
    step = (1000, 2.00, 1.50)
    if bounds.lower() == "kurucz95":
        bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]

    res, slopeEP, slopeRW, abundances = func(x0)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)
    if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
        return parameters, True

    try:
        os.remove('minimization_profile.dat')
    except OSError:
        pass

    N = 0
    while N < iteration:
        Nsub = 0
        cycle = [parameters[0]]
        while (abs(slopeEP) >= EPcriteria) and not fix_teff and Nsub < 15 and N < iteration:
            # For Teff
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**3
            if abs(step_i) < 1:
                step_i = s*1  # Minimum 1K
            elif abs(step_i) > 17625:
                step_i = s*17625  # Maximum 17625K (half the interval?)
            parameters[0] += step_i
            parameters[0] = check_bounds(parameters[0], bounds, 1)
            parameters[0] = int(parameters[0])
            if parameters[0] in cycle:
                break
            else:
                cycle.append(parameters[0])
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters, weights=weights)
            N += 1
            Nsub += 1
            Abdiff = np.diff(abundances)[0]
            save_iteration(parameters)
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

        Nsub = 0
        cycle = [parameters[3]]
        while (abs(slopeRW) >= RWcriteria) and not fix_vt and Nsub < 15 and N < iteration:
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
            parameters[3] = round(parameters[3], 2)
            if parameters[3] in cycle:
                break
            else:
                cycle.append(parameters[3])
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters, weights=weights)
            N += 1
            Nsub += 1
            Abdiff = np.diff(abundances)[0]
            save_iteration(parameters)
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

        Nsub = 0
        cycle = [parameters[1]]
        while (abs(Abdiff) >= ABdiffcriteria) and not fix_logg and Nsub < 15 and N < iteration:
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
            parameters[1] = round(parameters[1], 2)
            if parameters[1] in cycle:
                break
            else:
                cycle.append(parameters[1])
            print_format(parameters)
            res, slopeEP, slopeRW, abundances = func(parameters, weights=weights)
            N += 1
            Nsub += 1
            Abdiff = np.diff(abundances)[0]
            save_iteration(parameters)
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

        #       Input metal...   FeI abund.
        Nsub = 0
        while (parameters[2] != abundances[0]-7.47) and not fix_feh and Nsub < 5 and N < iteration:
            # For metalicity
            parameters[2] = abundances[0]-7.47
            parameters[2] = check_bounds(parameters[2], bounds, 5)
            parameters[2] = round(parameters[2], 2)
            res, slopeEP, slopeRW, abundances = func(parameters, weights=weights)
            N += 1
            Nsub += 1
            print_format(parameters)
            save_iteration(parameters)
        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

    print 'Stopped in %i iterations' % N
    converged = check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0])

    return parameters, converged


def _bump(x):
    """Bump to the values in the list, x"""
    for i, xi in enumerate(x):
        x[i] = np.random.normal(xi, abs(0.01*xi))
    return x


def minimize2(x0, func, bounds="kurucz95", weights='null',
             fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
             iteration=160, EPcrit=0.001, RWcrit=0.003, ABdiffcrit=0.01):
    """
    Sane minimization like a normal human being would do it.
    """
    global EPcriteria, RWcriteria, ABdiffcriteria
    global f_teff, f_logg, f_vt, f_feh
    EPcriteria, RWcriteria, ABdiffcriteria = EPcrit, RWcrit, ABdiffcrit
    f_teff, f_logg, f_feh, f_vt = fix_teff, fix_logg, fix_feh, fix_vt
    # Step size in Teff, logg, vt
    step = (700, 1.50, 0.50)
    if bounds.lower() == "kurucz95":
        bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]

    res, slopeEP, slopeRW, abundances = func(x0)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)
    if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
        return parameters, True
    try:
        os.remove('minimization_profile.dat')
    except OSError:
        pass

    all_params = [copy(parameters)]
    N = 0
    while N < iteration:
        # Step for Teff
        if (abs(slopeEP) >= EPcriteria) and not fix_teff:
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**3
            if abs(step_i) < 1:
                step_i = s*1  # Minimum 1K
            elif abs(step_i) > 17625:
                step_i = s*17625  # Maximum 17625K (half the interval?)
            parameters[0] += step_i
            parameters[0] = check_bounds(parameters[0], bounds, 1)
            parameters[0] = int(parameters[0])

        # Step for VT
        if (abs(slopeRW) >= RWcriteria) and not fix_vt:
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
            parameters[3] = round(parameters[3], 2)

        # Step for logg
        if (abs(Abdiff) >= ABdiffcrit) and not fix_logg:
            s = -np.sign(Abdiff)
            step_i = s * step[1]/abs(np.log(abs(Abdiff)+0.0005))**3
            if abs(step_i) < 0.01:
                step_i = s*0.01  # Minimum logg step
            elif abs(step_i) > 2.5:
                step_i = s*2.5  # Maximum logg step
            parameters[1] += step_i
            parameters[1] = check_bounds(parameters[1], bounds, 3)
            parameters[1] = round(parameters[1], 2)

        # Step for [Fe/H]
        if not fix_feh:
            parameters[2] = abundances[0]-7.47
            parameters[2] = check_bounds(parameters[2], bounds, 5)
            parameters[2] = round(parameters[2], 2)

        if parameters in all_params:
            parameters = _bump(parameters)
            parameters[0] = int(parameters[0])
            parameters[1] = round(parameters[1], 2)
            parameters[2] = round(parameters[2], 2)
            parameters[3] = round(parameters[3], 2)
        all_params.append(copy(parameters))

        res, slopeEP, slopeRW, abundances = func(parameters, weights=weights)
        Abdiff = np.diff(abundances)[0]
        N += 1
        print_format(parameters)
        save_iteration(parameters)

        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

    print 'Stopped in %i iterations' % N
    c = check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0])
    return parameters, c
