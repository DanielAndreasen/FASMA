#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from glob import glob
import os
import logging
# from model_interpolation import interpolator, save_model
from model_interpolation_new import interpolator, save_model


# Why a single leading underscore? Because this is an internal function.
# See pep8 for more information here:
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

    assert (t_rng[0] <= teff <= t_rng[1]), 'Teff out of range: %s to %s' % t_rng
    assert (l_rng[0] <= logg <= l_rng[1]), 'logg out of range: %s to %s' % l_rng
    assert (f_rng[0] <= feh <= f_rng[1]), '[Fe/H] out of range: %s to %s' % f_rng

    # Make the slice in [Fe/H]
    folders = glob('kurucz95/m*') + glob('kurucz95/p*')
    feh_grid = [folder.replace('kurucz95/', '') for folder in folders]
    feh_grid = [v.replace('m', '-') if v.startswith('m') else
                v.replace('p', '') for v in feh_grid]

    feh_grid = np.array(list(map(float, feh_grid))) / 10

    feh_m = feh_grid[abs(feh_grid - feh).argsort()[:2]]
    feh_m = [str(f).replace('.', '') for f in feh_m]
    paths = [f.replace('-', 'm') if f.startswith('-') else
             'p' + f for f in feh_m]

    # All the models in the [Fe/H]-space
    models = []
    for path in paths:
        models.extend(glob('kurucz95/%s/*.gz' % path))
    models = np.array(models)

    # This is a bit complicated way to get the temp. from the path of
    # all the models
    teff_m = [int(model.split('/')[-1].split('g')[0]) for model in models]
    diff_teff = abs(np.array(teff_m) - teff)
    idx_teff = []
    # Get the temperature closest and second closest to the teff selected. If
    # third closest is also needed, increace the loop by 1.
    for i in range(2):
        idx = np.where(diff_teff == min(diff_teff))[0]
        diff_teff[idx] = 99999
        idx_teff += list(idx)
    models = models[idx_teff]

    logg_m = [model.replace(path, '').split('g')[1].split('.')[0] for model in models]
    logg_m = np.array(map(float, logg_m)) / 10
    diff_logg = abs(logg_m - logg)
    idx_logg = []
    for i in range(2):
        idx = np.where(diff_logg == min(diff_logg))[0]
        diff_logg[idx] = 99
        idx_logg += list(idx)

    nn_feh = sorted(tuple(float(f[0:-1]+'.'+f[-1]) for f in feh_m))
    nn_logg = sorted(tuple(set(logg_m[idx_logg])))
    nn_teff = sorted(tuple(set(np.array(teff_m)[idx_teff])))

    models = sorted(models[idx_logg])
    return models, nn_teff, nn_logg, nn_feh,


# TODO: This function will be merged with the one beneath
def _update_par(infile='batch.par', out='batch.par'):
    """Update the parameter file

    :infile: The input file
    :out: The output file
    """

    '''
    Runs MOOG with the given input parameters and returns a numpy array of
    the outputted smooth spectrum.

    Inputs
    -----

    atmosphere_model    :   Location of your model atmosphere file
    line_list           :   Location of your line list

    Additional keyword arguments
    ----------------------------

    These additional keyword arguments allow the user to have full control
    over what is put into the MOOG input file. The default values are:

    terminal        'x11'
    atmosphere      1
    molecules       2
    trudamp         1
    lines           1
    flux/int        1
    damping         0
    units           0
    iraf            0
    plot            2
    obspectrum      1       Unless obspectrum is provided to the function.
    opacit          0
    freeform        0
    strong          0       Unless a strong lines list is provided.
    plotpars        1       0.75 Gaussian smoothing by default. Show full
                            synthesized spectral range with y:[0, 1.2]
    histogram       0
    synlimits               Defaults to the wavelength range provided and
                            the given wavelength step size, and the delta
                            defaults to the wavelength step size.

    Outputs
    -------

    A numpy two-dimensional spectrum which contains the wavelength in the
    entire first column, and normalised smoothed flux in the second column
    '''


def run(atmosphere_model='out.atm', line_list='linelist.moog', **kwargs):

    # Path checks for input files
    if not os.path.exists(atmosphere_model):
        raise IOError('Atmosphere model file "%s" could not be found.' %
                      (atmosphere_model))

    if not os.path.exists(line_list):
        raise IOError('Line list file "%s" could not be found.' %
                      (line_list))

    default_kwargs = {
        'atmosphere': 1,
        'molecules':  1,
        'trudamp':    1,  # Sure, why not? It's a black art anyway!
        'lines':      1,
        'terminal':   'x11',
        'flux/int':   0,
        'damping':    2,
        'units':      0,
        'iraf':       0,
        'plot':       0,
        'obspectrum': 0,
        'opacit':     0,
        'freeform':   0,
        'strong':     0,
        }

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.iteritems():
        if key not in kwargs.keys():
            kwargs[key] = value

    # Generate a MOOG-compatible run file
    moog_filename = 'batch.par'

    if os.path.exists(moog_filename):
        logging.warn('Temporary MOOG file already exists (%s), over-writing...'
                     % moog_filename)

    moog_contents = "abfind\n"\
                    "terminal       %s\n"\
                    "model_in       '%s'\n"\
                    "summary_out    '%s'\n"\
                    "standard_out   '%s'\n"\
                    "lines_in       '%s'\n"  % (kwargs['terminal'], atmosphere_model, 'summary.out',
                           'result.out', line_list)

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,plot,opacit,freeform,obspectrum,histogram,'\
               'synlimits'.split(',')
    if kwargs.has_key('plotpars'):
        if kwargs['plotpars'] != 0:
            logging.warn('Plotting require SM')
            settings.append('plotpars')

    for setting in settings:
        if kwargs.has_key(setting):
            moog_contents += "%s %s\n" % (setting + ' ' * (14 - len(setting)),
                                        kwargs[setting])

    with open(moog_filename, 'w') as moog:
        moog.writelines(moog_contents)


def _transform_micro(teff, logg, feh):
    """
    This is the python equivalent to transform.f
    """
    # puts('Please type in the new ' + colored.yellow('micro turbulence'))
    v_micro = raw_input('> ')
    try:
        v_mico = float(v_micro)
    except ValueError:
        # puts(colored.red('Please provide a number...'))
        v_micro = _transform_micro(teff, logg, feh)
    return v_micro


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    teff, logg, feh = 5001, 4.01, 0.01
    models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=feh, type='kurucz95')

    # TODO: First thing to vary should be Fe/H, then logg, then Teff (see
    # around L450 in atlas9.py by jobovy)
    # for m in models:
        # print(m)

    # raise SystemExit('Exiting...')
    m_all, m_out, _ = interpolator(models, teff=(teff, nt), logg=(logg, nl),
                                   feh=(feh, nf))

    # m_all, m_out, _ = interpolator(models, teff=(teff, sorted(nt)[1:3]), logg=(logg, nl),
    #                                feh=(feh, nf))
    # save_model(m_out, vt=2.4)
