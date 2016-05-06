#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import pandas as pd
import argparse
from utils import _update_par as updateBatch
from utils import _run_moog as runMoog
from utils import Readmoog, GetModels
from interpolation import interpolator, save_model


def recalSingleLine(line, params=None, version=2014, maxiter=40):
    '''Recalibrate a single line and return the new loggf

    Inputs
    ------
    line : list
      The line containing (wavelength, element, EP, loggf, EW) in that order
    params : list/tuple
      The parameters (Teff, logg, [Fe/H], vt)
    version : int
      The version of MOOG

    Output
    ------
    loggf : float
      The new recalibrated loggf
    '''

    def moogAbund(loggf):
        line[3] = loggf
        np.savetxt('temporary.moog', line[:, np.newaxis].T, fmt=fmt, header=header)
        runMoog()
        m = Readmoog(params=params, version=version)
        _, abund = m.elements()
        return abund[0] - 7.47

    fmt = ('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f')
    header = 'Wavelength     ele    EP     loggf        EW'
    np.savetxt('temporary.moog', line[:, np.newaxis].T, fmt=fmt, header=header)
    loggf_old = line[3]
    a, b = -10, 10  # extreme values of loggf
    c = (a+b)/2
    for _ in range(maxiter):
        if c == 0:  # Don't evaluate at loggf = 0
            c += (abs(a) + abs(b)) / 10
        fa = moogAbund(a)
        fc = moogAbund(c)
        if fc == 0:
            return c
        elif fa*fc < 0:
            b = c
        else:
            a = c
        c = (a+b)/2
    return c


def _parser():
    parser = argparse.ArgumentParser(description='Recalibrate the loggf for a set of settings')
    parser.add_argument('input', help='Input linelist from rawLinelist folder')
    parser.add_argument('output', help='Name of output file (saved in rawLinelist)')
    parser.add_argument('-m', '--model', help='Model atmosphere', default='kurucz95', choices=['kurucz95', 'apogee_kurucz', 'marcs'])
    parser.add_argument('-v', '--moogversion', help='MOOG version', default=2014)
    parser.add_argument('-d', '--damping', help='Damping to be used in MOOG', default=1, choices=map(str, [1, 2]))
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = _parser()

    fname = args.input
    fout = 'rawLinelist/%s' % args.output
    lines = pd.read_csv(fname, skiprows=2, delimiter=r'\s+', names=['WL', 'num', 'EP', 'loggf', 'ele', 'EW'])

    params = [5777, 4.44, 0.00, 1.00]
    teff, logg, feh, _ = params
    grid = GetModels(teff, logg, feh, args.model)
    models, nt, nl, nf = grid.getmodels()
    model = interpolator(models, teff=(teff, nt), logg=(logg, nl), feh=(feh, nf))
    save_model(model, params)

    cols = ['WL', 'num', 'EP', 'loggf', 'EW']
    fmt = ('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f')
    header = 'Wavelength     ele    EP     loggf        EW'
    x = lines[cols].values[0][:, np.newaxis].T
    np.savetxt('temporary.moog', x, fmt=fmt, header=header)

    options = {'driver': 'abfind',
               'damping': args.damping}
    updateBatch(line_list='temporary.moog', **options)

    newloggf = np.zeros(lines.shape[0])
    for i, line in enumerate(lines[cols].values):
        print 'Wavelength: %.3f' % line[0]
        print 'Old loggf: %.3f' % line[3]
        zz = recalSingleLine(line, params=params, version=args.moogversion)
        zz = np.log10(zz) if zz > 0 else zz
        newloggf[i] = zz
        print 'New loggf: %.3f\n' % newloggf[i]

    lines['newloggf'] = pd.Series(newloggf)
    X = lines[['WL', 'num', 'EP', 'newloggf', 'EW']]
    print 'Saving results to: %s' % fout
    np.savetxt(fout, X, fmt=fmt, header=header)
