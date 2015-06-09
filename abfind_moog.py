#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

# Path to the interpolation code
path = '/home/daniel/Software/SPECPAR/interpol_models/'


def _create_model(teff, logg, feh, vmicro):
    """Create the atmosphere models"""
    # Interpolate over models
    intermod_par = tuple(map(str, [teff, logg, feh]) + [path])
    cmd = 'echo %s %s %s | %sintermod.e > tmp' % (intermod_par)
    os.system(cmd)

    # Transform to desired micro turbulence
    cmd = 'echo %s | %stransform.e > tmp' % (str(vmicro), path)
    os.system(cmd)

    # Clean after
    os.system('rm -f mod? fort* tmp')


def _update_batch(linelist):
    """Update the batch.par file for MOOG"""
    with open('batch.par', 'r') as lines:
        lines_tmp = []
        for line in lines:
            if 'lines_in' in line:
                line = line.split('      ')
                line[1] = "'" + linelist + "'\n"
                line = '      '.join(line)
            lines_tmp.append(line)

    with open('batch.par', 'w') as file:
        file.write(''.join(lines_tmp))


def _data_structure(out='summary.log'):
    """Black magic here"""
    structure = []
    with open(out, 'r') as lines:
        s = 0
        for line in lines:
            line = filter(None, line.strip().split())
            try:
                line = map(float, line)
            except ValueError, e:
                if s > 0:
                    structure.append(s)
                line = []
                s = 0
            if len(line) > 1:
                s += 1
    return structure


def read_output(out='summary.log'):
    """
    Read the summary file and extract the elements and return a numpy structure
    of the data.

    nelements is the number of elements, e.g. 2 if analysing FeI and FeII.
    """

    structure = _data_structure(out=out)
    nelements = len(structure)

    element = [[]] * nelements
    for i, struct in enumerate(structure):
        element[i] = np.zeros((struct, 7))

    i = -1  # counter for which element we are saving
    with open(out, 'r') as lines:
        s = -1
        for line in lines:
            s += 1
            if line.startswith('wavelength'):
                i += 1

            line = filter(None, line.strip().split())
            try:
                line = map(float, line)
            except ValueError, e:
                line = []
                s = -1

            if len(line) > 1:
                element[i][s] = line

    return element


def plot_data(data):
    EP = data[:, 1]
    logRW = data[:, 4]
    abund = data[:, 5]
    m = np.mean(abund)
    s = np.std(abund)

    z1 = np.polyfit(EP, abund, 1)
    p1 = np.poly1d(z1)

    z2 = np.polyfit(logRW, abund, 1)
    p2 = np.poly1d(z2)

    plt.figure(figsize=(8, 9))
    plt.subplot(211)
    plt.plot(EP, abund, '.w', label='Temperature')
    plt.plot(EP, abund, '.k')
    plt.hlines(m, min(EP), max(EP), colors='b', linestyles='--',
               linewidth=3)
    plt.hlines(m - s, min(EP), max(EP), colors='b', linestyles='--',
               linewidth=3)
    plt.hlines(m + s, min(EP), max(EP), colors='b', linestyles='--',
               linewidth=3)
    plt.plot(EP, p1(EP), '-r')
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Excitation potential: $\chi$ [eV]')
    plt.ylabel('Abundance')

    plt.subplot(212)
    plt.plot(logRW, abund, '.w', label='Surface gravity')
    plt.plot(logRW, abund, '.k')
    plt.hlines(m, min(logRW), max(logRW), colors='b', linestyles='--',
               linewidth=3)
    plt.hlines(m - s, min(logRW), max(logRW), colors='b', linestyles='--',
               linewidth=3)
    plt.hlines(m + s, min(logRW), max(logRW), colors='b', linestyles='--',
               linewidth=3)
    plt.plot(logRW, p2(logRW), '-r')
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Reduced EW: $\log RW$')
    plt.ylabel('Abundance')
    return z1[0], z2[0]


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Force fit abundances with MOOG',
                                epilog='Happy spectroscopying :)')

    p.add_argument('teff', type=int, help='The effective temperature')
    p.add_argument('logg', type=float, help='The surface gravity')
    p.add_argument('feh', type=float, help='The metallicity')
    p.add_argument('vmicro', type=float, help='The micro turbulence')
    p.add_argument('-l', '--linelist', help='The linelist to be used')
    p.add_argument('-o', '--output', help='The output file with abundances',
                   default='summary.log')
    p.add_argument('-p', '--plot', help='Enable plotting', default=True,
                   type=bool)

    args = p.parse_args()

    # Interpolate and transform models
    _create_model(args.teff, args.logg, args.feh, args.vmicro)

    # Update the batch file if a new linelist is used
    if args.linelist:
        _update_batch(args.linelist)

    # Run moog
    os.system('MOOGSILENT > zzz')
    os.system('rm -f zzz')

    # Prepare the data
    data = read_output(args.output)

    # plt.ion()
    c1, c2 = plot_data(data[0])
    if c1 > 0:
        print('Raise temperature')
    elif c1 < 0:
        print('Lower temperature')
    else:
        print('Temperature is found %i' % args.temperature)
    raw_input('Press RETURN to exit > ')

    if args.plot:
        plt.ion()
        for d in data:
            plot_data(d)
        raw_input('Press RETURN to exit > ')
