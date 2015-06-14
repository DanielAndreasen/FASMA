#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import os
import argparse
import seaborn as sns
sns.set_style('dark')
sns.set_context('talk')
import matplotlib.pyplot as plt

from pymoog import _get_model
from model_interpolation import interpolator
from model_interpolation import save_model
from moog_minimization import minimize
from runmoog import fun_moog, fun_moog_fortran



def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Plot fits file for ARES. Be'
                                     ' careful with large files')


    parser.add_argument('linelist', help='Input line list')
    parser.add_argument('-p', '--parfile',
            help='The parameter file (default: batch.par)',
            default='batch.par')
    parser.add_argument('-m', '--model',
            help='Model atmosphere',
            default='Kurucz95',
            choices=['Kurucz95', 'Kn', 'Marcs', 'PHOENIX'])
    parser.add_argument('-i', '--initial',
            help='Initial conditions (Teff, logg, [Fe/H], vt)',
            nargs='+',
            type=float,
            default=False)
    parser.add_argument('-f', '--fix',
            help='Parameters to fix',
            nargs='+',
            type=int,
            default=[0, 0, 0, 0])
    parser.add_argument('-pl', '--plot',
            help='Plot the slopes',
            default=False)
    parser.add_argument('-ol', '--outliers',
            help='Remove n*sigma outliers',
            type=float,
            default=False)
    parser.add_argument('-spt', '--spectralType',
            help='Input spectral type (e.g. F4V) and get initial parameters',
            default=False)
    parser.add_argument('-v', '--verbose',
            help='Print information to the screen along the way',
            default=False)
    args = parser.parse_args()
    return args


def moogme(linelist, parfile='batch.par', model='Kurucz95',
           initial=False, fix_params=(0, 0, 0, 0),
           plot=False, outlier=False, spt=False):
    """
    Some doc
    """

    # TODO: Someone put the right temperatures and logg below. Please
    # keep the structure even though the lines are longer than 80 chars.
    # Temperatures for V from http://www.uni.edu/morgans/astro/course/Notes/section2/spectraltemps.html
    spectralType_T = {                                                                'O5': 54000, 'O6': 45000, 'O7': 43300, 'O8': 40600, 'O9': 37800,
                      'B0': 29200, 'B1': 23000, 'B2': 21000, 'B3': 17600,             'B5': 15200, 'B6': 14300, 'B7': 13500, 'B8': 12300, 'B9': 11400,
                      'A0': 9600,  'A1': 9330,  'A2': 9040,  'A3': 8750,  'A4': 8480, 'A5': 8310,               'A7': 7920,
                      'F0': 7350,               'F2': 7050,  'F3': 6850,              'F5': 6700,  'F6': 6550,  'F7': 6400,  'F8': 6300,
                      'G0': 6050,  'G1': 5930,  'G2': 5800,                           'G5': 5660,                            'G8': 5440,
                      'K0': 5240,  'K1': 5110,  'K2': 4960,  'K3': 4800,  'K4': 4600, 'K5': 4400,               'K7': 4000,
                      'M0': 3750,  'M1': 3700,  'M2': 3600,  'M3': 3500,  'M4': 3400, 'M5': 3200,  'M6': 3100,  'M7': 2900,  'M8': 2700}
    # Daniel's late-night loggs. Maybe this should be changed. It is 02:12 AM saturday, and I'm tired...
    spectralType_g = {'I': 0, 'II': 1, 'III': 2, 'IV': 3, 'V': 4.5}

    # Setting the initial parameters
    if spt and not initial:
        Teff = spectralType_T[spt[0:2]]
        logg = spectralType_g[spt[2::]]
        initial = (Teff, logg, 0.00, 1.00)
    elif not initial and not spt:
        # Set initial parameters to solar value
        initial = (5777, 4.44, 0.00, 1.00)

    print(initial)

    # Preparing the batch file for MOOGSILENT
    if not os.path.isfile(parfile):
        # TODO: Maybe we can create one. The function is in pymoog.py for that
        raise IOError('%s does not exist' % parfile)
    if parfile != 'batch.par':
        os.system('cp %s batch.par' % parfile)
        rm_batch = True


    # Setting the models to use
    if model != 'Kurucz95':
        raise NotImplementedError('Your request for type: %s is not available' % model)

    # Setting which parameters to fix
    fix_teff = True if fix_params[0] else False
    fix_logg = True if fix_params[1] else False
    fix_feh  = True if fix_params[2] else False
    fix_vt   = True if fix_params[3] else False



    # Get the initial grid models
    # models, nt, nl, nf = _get_model(teff=initial[0], logg=initial[1], feh=initial[2])
    # inter_model = interpolator(models,
    #                            teff=(initial[0], nt),
    #                            logg=(initial[1], nl),
    #                            feh=(initial[2], nf))
    # save_model(inter_model, params=initial)

    # parameters = minimize(initial, fun_moog_fortran, bounds=model,
    #                       fix_teff=fix_teff, fix_logg=fix_logg,
    #                       fix_feh=fix_feh, fix_vt=fix_vt)

    # print('Congratulation, you have won! Your final parameters are\n' + parameters)
    return 0
    # return parameters


if __name__ == '__main__':
    args = _parser()
    # print(args)
    parameters = moogme(args.linelist, args.parfile, args.model,
                        args.initial, args.fix, args.plot,
                        args.outliers, args.spectralType)
