#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from glob import glob
# from clint.textui import puts, colored, indent


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

    # Let's stay in range of the models available and do some checks
    if type.lower() == 'kurucz95':
        t_rng = 3500, 35000
        l_rng = 2.5, 5.0
        f_rng = -3.0, 1.0
    elif type.lower() == 'marcz':  # How to spell this?
        raise NotImplementedError('Patience is the key. Wait a bit more for %s\
                                   models to be implemented.' % type)
    else:
        raise NameError('Could not find %s models' % type)
    assert (t_rng[0] <= teff and teff <= t_rng[1]), 'Teff out of range: %s to %s' % t_rng
    assert (l_rng[0] <= logg and logg <= l_rng[1]), 'logg out of range: %s to %s' % l_rng
    assert (f_rng[0] <= feh and feh <= f_rng[1]), '[Fe/H] out of range: %s to %s' % f_rng

    # Make the slice in [Fe/H]
    folders = glob('kurucz95/m*') + glob('kurucz95/p*')
    feh_grid = [folder.replace('kurucz95/', '') for folder in folders]
    feh_grid = [v.replace('m', '-') if v.startswith('m') else
                v.replace('p', '') for v in feh_grid]
    feh_grid = np.sort(np.array(map(float, feh_grid)) / 10)

    feh_m = []
    for _ in range(2):
        i = np.argmin(abs(feh_grid - feh))
        feh_m.append(str(feh_grid[i]).replace('.', ''))
        feh_grid[i] = 9

    paths = [f.replace('-', 'm') if f.startswith('-') else
             'p' + f for f in feh_m]

    # All the models in the [Fe/H]-space
    models = []
    for path in paths:
        models.extend(glob('kurucz95/%s/*.gz' % path))
    models = np.array(models)

    # This is a bit complicated way to get the temp. from the path of all the
    # models
    teff_m = [int(model.split('/')[-1].split('g')[0]) for model in models]
    diff_teff = abs(np.array(teff_m) - teff)
    idx_teff = []
    # Get the temperature closest and second closest to the teff selected. If
    # third closest is also needed, increace the loop by 1.
    for i in range(4):
        idx = np.where(diff_teff == min(diff_teff))[0]
        diff_teff[idx] = 99999
        idx_teff += list(idx)
    models = models[idx_teff]

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
    # puts('Please type in the new ' + colored.yellow('micro turbulence'))
    v_micro = raw_input('> ')
    try:
        v_mico = float(v_micro)
    except ValueError, e:
        # puts(colored.red('Please provide a number...'))
        v_micro = _transform_micro(teff, logg, feh)
    return v_micro


def _interpolate_atm():
    pass


def minimize(teff, logg, feh):
    pass


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    # from sys import argv
    print np.sort(_get_model(teff=5777, logg=4.44, feh=0.00, type='kurucz95'))
    # _transform_micro(3750, 4.40, -0.5)
