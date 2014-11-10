#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from glob import glob
from clint.textui import puts, colored, indent


# Why a single leading underscore? Because this is an internal function. See
# pep8 for more information here:
# http://legacy.python.org/dev/peps/pep-0008/#naming-conventions
def _run_moog():
    os.system('MOOGSILENT > tmp')
    return 0


def _get_model(teff, logg, feh, type='kurucz95'):
    """
    Get the path to the model with the given effective temperature, logg and
    metallicity
    """

    # Let's keep in range of the models available
    assert (35000 >= teff and teff >= 3500), 'Temperature out of range'
    assert (5.0 >= logg and logg >= 2.5), 'Surface gravity out of range'

    if feh < 0:
        sign = '/m'
        # Remove the '-' and the '.'
        feh = str(feh).replace('-', '').replace('.', '')
    else:
        sign = '/p'
        feh = str(feh).replace('.', '')
    folder = sign + str(feh) + '/'
    path = type + folder
    models = np.array(glob(path + '*.gz'))
    # This is a bit complicated way to get the temp. from the path of all the
    # models
    teff_m = [int(model.replace(path, '').split('g')[0]) for model in models]
    diff_teff = abs(np.array(teff_m) - teff)
    idx_teff = []
    # Get the temperature closest and second closest to the teff selected. If
    # third closest is also needed, increace the loop by 1.
    for i in range(2):
        idx = np.where(diff_teff == min(diff_teff))[0]
        diff_teff[idx] = 99999
        idx_teff += list(idx)
    models = models[idx_teff]
    # print models


    logg_m = [model.replace(path, '').split('g')[1].split('.')[0] for model in models]
    logg_m = np.array([float(li)/10 for li in logg_m])
    diff_logg = abs(logg_m - logg)
    idx_logg = []
    for i in range(2):
        idx = np.where(diff_logg == min(diff_logg))[0]
        diff_logg[idx] = 99
        idx_logg += list(idx)

    models = models[idx_logg]
    return models



def _transform_micro(teff, logg, feh):
    """
    This is the python equivalent to transform.f
    """
    puts('Please type in the new ' + colored.yellow('micro turbulence'))
    v_micro = raw_input('> ')
    try:
        v_mico = float(v_micro)
    except ValueError, e:
        puts(colored.red('Please provide a number...'))
        v_micro = _transform_micro(teff, logg, feh)
    return v_micro


def _interpolate_atm():
    pass


def minimize(teff, logg, feh):
    pass


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    # from sys import argv
    print _get_model(6777, 4.4, 0.0)
    # _transform_micro(3750, 4.40, -0.5)


