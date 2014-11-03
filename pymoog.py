#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from glob import glob


# Why a single leading underscore? Because this is an internal function. See
# pep8 for more information here:
# http://legacy.python.org/dev/peps/pep-0008/#naming-conventions
def _run_moog():
    os.system('MOOGSILENT > tmp')
    return 0


def _get_model(teff, logg, feh, type='kurucz95'):
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
    models = glob(path + '*.gz')
    for model in models:
        model = model.replace(path, '')
        model_m = model.split('g')
        teff_m = int(model_m[0])
        if teff_m == teff:
            logg_m = model_m[1].split(sign[1])[0]
            logg_m = float(logg_m.strip('.'))/10
            if logg_m == logg:
                return model


def _interpolate_atm():
    return 0


def minimize(teff, logg, feh):
    return 0


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    # from sys import argv
    print _get_model(3750, 3.0, 0.1)


