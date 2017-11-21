#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import pandas as pd
import argparse
from utils import _update_par as updateBatch
from utils import _run_moog as runMoog
from utils import Readmoog
from interpolation import interpolator
import os


def solar_abundance(atom):
    '''Give atomic number and return solar abundance from Asplund et al. 2009

    Input
    -----
    atom : int
      The atomic number

    Output
    ------
    abundance : float
      The solar abundance of the atom
    '''
    if not isinstance(atom, int):
        raise ValueError('Atomic number need to be an integer')
    solar = [12.00, 10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,
             6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,
             3.15, 4.95, 3.93, 5.64, 5.43, 7.47, 4.99, 6.22, 4.19, 4.56,
             3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,
             1.46, 1.88, -5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,
             1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,
             -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,
             0.10, 0.85, -0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,
             0.90, 1.75, 0.65, -5.00, -5.00, -5.00, -5.00, -5.00, -5.00,
             0.02, -5.00, -0.54, -5.00, -5.00, -5.00]
    return solar[atom-1]


def recalSingleLine(line, params=None, version=2014, maxiter=40, driver='abfind'):
    '''Recalibrate a single line and return the new loggf

    Inputs
    ------
    line : list
      The line containing (wavelength, element, EP, loggf, EW) in that order
    params : list/tuple
      The parameters (Teff, logg, [Fe/H], vt)
    version : int
      The version of MOOG
    driver : str
      The MOOG driver to use (abfind or ewfind)

    Output
    ------
    loggf : float
      The new recalibrated loggf
    '''

    def moogAbund(loggf, ewdriver=False):
        line[3] = loggf
        np.savetxt('temporary.moog', line[:, np.newaxis].T, fmt=fmt, header=header)
        runMoog()
        if ewdriver:
            d = np.loadtxt('summary.out', skiprows=5, usecols=(6,))
            out = d-line[4]
        else:
            m = Readmoog(params=params, version=version)
            _, abund = m.elements()
            solar = solar_abundance(int(line[1]))
            out = round(abund[0] - solar, 3)
        return out

    ewdriver = True if driver == 'ewfind' else False

    fmt = ('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f')
    header = 'Wavelength     ele    EP     loggf        EW'
    np.savetxt('temporary.moog', line[:, np.newaxis].T, fmt=fmt, header=header)
    loggf_old = line[3]
    a, b = loggf_old-5, loggf_old+5  # extreme values of loggf
    c = (a+b)/2
    for _ in range(maxiter):
        if c == 0:  # Don't evaluate at loggf = 0
            c += (abs(a) + abs(b)) / 10
        fa = moogAbund(a, ewdriver=ewdriver)
        fc = moogAbund(c, ewdriver=ewdriver)
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
    parser.add_argument('-dr', '--driver', help='Which driver to use', default='abfind', choices=['abfind', 'ewfind'])
    parser.add_argument('-p', '--parameters', help='Atmospheric parameters, Teff, logg, [Fe/H], vt', nargs='+', default=None)
    args = parser.parse_args()
    return args


if __name__ == '__main__':

    args = _parser()

    fname = args.input
    fout1 = 'rawLinelist/%s' % args.output
    fout2 = 'linelist/%s' % args.output.replace('.ares', '.moog')
    lines = pd.read_csv(fname, skiprows=2, delimiter=r'\s+', names=['WL', 'num', 'EP', 'loggf', 'ele', 'EW'])

    if args.parameters is None:
        params = [5777, 4.44, 0.00, 1.00]
    else:
        params = map(float, args.parameters)
        params[0] = int(params[0])

    interpolator(params=params, atmtype=args.model, save=True)

    cols = ['WL', 'num', 'EP', 'loggf', 'EW']
    fmt2 = ('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f')
    header1 = 'WL         num       E.P.     loggf         ele     EWsun\n'
    header1 += '-------    ----      ----     ------        ----    -----'
    header2 = 'Wavelength     ele       EP      loggf                          EW'
    x = lines[cols].values[0][:, np.newaxis].T
    np.savetxt('temporary.moog', x, fmt=fmt2, header=header2)

    options = {'driver': args.driver,
               'damping': args.damping}
    updateBatch(line_list='temporary.moog', **options)

    newloggf = np.zeros(lines.shape[0])
    for i, line in enumerate(lines[cols].values):
        print 'Wavelength: %.3f' % line[0]
        print 'Old loggf: %.3f' % line[3]
        zz = recalSingleLine(line, params=params, version=args.moogversion, driver=args.driver)
        zz = np.log10(zz) if zz > 0 else zz
        newloggf[i] = zz
        print 'New loggf: %.3f\n' % newloggf[i]

    lines['newloggf'] = pd.Series(newloggf)
    X = lines[['WL', 'num', 'EP', 'newloggf', 'ele', 'EW']]
    fmt1 = ('%7.2f', '%7.1f', '%9.2f', '%10.3f', '%10s', '%9.1f')
    print 'Saving results to: %s' % fout1
    np.savetxt(fout1, X, fmt=fmt1, header=header1, comments='')

    X = lines[['WL', 'num', 'EP', 'newloggf', 'EW']]
    print 'Saving results to: %s' % fout2
    np.savetxt(fout2, X, fmt=fmt2, header=header2)
    os.remove('temporary.moog')
