#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from model_interpolation import interpolator, save_model

K95 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                32000, 33000, 34000, 35000, 3500, 36000, 37000, 38000, 39000),
       'feh': (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
               0.1, 0.2, 0.3, 0.5, 1.0),
       'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}


def _get_model(teff, logg, feh, type='kurucz95'):
    """
    Find the names of the closest grid points for a given effective
    temperature, surface gravity, and iron abundance (proxy for metallicity).

    Inputs
    -----

    teff  :   The effective temperature(K) for the model atmosphere
    logg  :   The surface gravity (logarithmic in cgs) for the model atmosphere
    feh   :   The metallicity for the model atmosphere
    type  :   The type of atmosphere models to use. Currently only Kurucz from
              '95 is supported.

    output
    ------
    models      : List with path to 8 models the two closest in each parameter
                  space (2x2x2)
    teff_model  : The two closest effective temperatures in the grid
    logg_model  : The two closest surface gravities in the grid
    feh_model   : The two closest metallicities in the grid

    The last three return values are used for the interpolation to do some
    mapping. If only the paths to the models are needed, do not pay attention
    to them.

    """

    if type == 'kurucz95':
        grid = K95
    else:
        raise NotImplementedError('You request for type: %s is not available' %
                                  type)
    if (teff < grid['teff'][0]) or (teff > grid['teff'][-1]):
        raise ValueError('Teff out of bounds: %s' % teff)
    if (logg < grid['logg'][0]) or (logg > grid['logg'][-1]):
        raise ValueError('logg out of bounds: %s' % logg)
    if (feh < grid['feh'][0]) or (feh > grid['feh'][-1]):
        raise ValueError('[Fe/H] out of bounds: %s' % feh)

    if teff in grid['teff']:
        teff_model = [teff, teff]
    else:
        for i, Ti in enumerate(grid['teff']):
            if Ti - teff > 0:
                break
        teff_model = [grid['teff'][i-1], grid['teff'][i]]

    if logg in grid['logg']:
        logg_model = [logg, logg]
    else:
        for i, li in enumerate(grid['logg']):
            if li - logg > 0:
                break
        logg_model = [grid['logg'][i-1], grid['logg'][i]]

    if feh in grid['feh']:
        feh_model = [feh, feh]
    else:
        for i, fi in enumerate(grid['feh']):
            if fi - feh > 0:
                break
        feh_model = [grid['feh'][i-1], grid['feh'][i]]

    name = lambda t, g, s, f: 'kurucz95/%s%s/%ig%i.%s%s.gz' % (s, f, t, g*10, s, f)
    models = []
    for teff_m in teff_model:
        for logg_m in logg_model:
            for feh_m in feh_model:
                if feh_m >= 0:
                    feh_m = str(feh_m).replace('.', '')
                    fout = name(teff_m, logg_m, 'p', feh_m)
                else:
                    feh_m = str(feh_m).replace('.', '').replace('-', '')
                    fout = name(teff_m, logg_m, 'm', feh_m)
                models.append(fout)

    return models, teff_model, logg_model, feh_model


def _update_par(atmosphere_model='out.atm', line_list='linelist.moog',
                infile='batch.par', **kwargs):
    """Runs MOOG with the given input parameters and returns a numpy array of
    the outputted smooth spectrum.

    Inputs
    -----

    atmosphere_model    :   Location of your model atmosphere file
    line_list           :   Location of your line list
    infile              :   Name of the parameter file

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

    And updated parameter file
    """

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

    moog_contents = "abfind\n"\
                    "terminal       %s\n"\
                    "model_in       '%s'\n"\
                    "summary_out    '%s'\n"\
                    "standard_out   '%s'\n"\
                    "lines_in       '%s'\n" % (kwargs['terminal'],
                                               atmosphere_model, 'summary.out',
                                               'result.out', line_list)

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,plot,opacit,freeform,obspectrum,histogram,'\
               'synlimits'.split(',')
    if 'plotpars' in kwargs:
        if kwargs['plotpars'] != 0:
            settings.append('plotpars')

    for setting in settings:
        if setting in kwargs:
            moog_contents += "%s %s\n" % (setting + ' ' * (14 - len(setting)),
                                          kwargs[setting])

    with open(moog_filename, 'w') as moog:
        moog.writelines(moog_contents)


def main(teff, logg, feh):
    models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=feh)
    n = interpolator(models, teff=(teff, nt), logg=(logg, nl), feh=(feh, nf))
    save_model(n, params=(teff, logg, feh, 2.4))


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    teff, logg, feh = 4252, 3.52, -0.02
    main(teff, logg, feh)
