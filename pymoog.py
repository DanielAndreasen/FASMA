#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from glob import glob
import os
import logging
from model_interpolation import interpolator, save_model

K95_teff = (3750,4000,4250,4500,4750,5000,5250,5500,5750,6000,
        6250,6500,6750,7000,7250,7500,7750,8000,8250,8500,8750,9000,9250,9500,
        9750,10000,10250,10500,10750,11000,11250,11500,11750,12000,12250,12500,
        12750,13000,14000,15000,16000,17000,18000,19000,20000,21000,22000,23000,
        24000,25000,26000,27000,28000,29000,30000,31000,32000,33000,34000,
        35000,3500,36000,37000,38000,39000)
K95_feh = (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2,
        0.3, 0.5, 1.0)
K95_logg = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)


def _get_model(teff, logg, feh, type='kurucz95'):
    if (teff < K95_teff[0]) or (teff > K95_teff[-1]):
        raise ValueError('Teff out of bounds: %s' % teff)
    if (logg < K95_logg[0]) or (logg > K95_logg[-1]):
        raise ValueError('logg out of bounds: %s' % logg)
    if (feh < K95_feh[0]) or (feh > K95_feh[-1]):
        raise ValueError('[Fe/H] out of bounds: %s' % feh)

    if teff in K95_teff:
        teff_model = [teff, teff]
    else:
        for i, K95T in enumerate(K95_teff):
            if K95T - teff > 0:
                break
        teff_model = [K95_teff[i-1], K95_teff[i]]

    if logg in K95_logg:
        logg_model = [logg, logg]
    else:
        for i, K95L in enumerate(K95_logg):
            if K95L - logg > 0:
                break
        logg_model = [K95_logg[i-1], K95_logg[i]]

    if feh in K95_feh:
        feh_model = [feh, feh]
    else:
        for i, K95F in enumerate(K95_feh):
            if K95F - feh > 0:
                break
        feh_model = [K95_feh[i-1], K95_feh[i]]

    name = lambda t, g, s, f: 'kurucz95/%s%s/%sg%i.%s%s.gz' % (s, f, t, g*10, s, f)
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

    # print models
    return models, teff_model, logg_model, feh_model


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


if __name__ == '__main__':
    # This is only for testing and should be removed later on...
    teff, logg, feh = 4250, 3.50, -0.02
    models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=feh)
    n = interpolator(models, teff=(teff, nt), logg=(logg, nl), feh=(feh, nf))
    save_model(n, params=(teff, logg, feh, 2.4))
