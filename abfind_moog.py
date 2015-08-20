#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import seaborn as sns
sns.set_style('dark')
sns.set_context('notebook', font_scale=1.5)
c = sns.color_palette()

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

    with open('batch.par', 'w') as lines:
        lines.write(''.join(lines_tmp))


def _data_structure(out='summary.out'):
    """Black magic here"""
    structure = []
    with open(out, 'r') as lines:
        s = 0
        for line in lines:
            line = filter(None, line.strip().split())
            try:
                line = map(float, line)
            except ValueError:
                if s > 0:
                    structure.append(s)
                line = []
                s = 0
            if len(line) > 1:
                s += 1
    return structure


def read_output(out='summary.out'):
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
            except ValueError:
                line = []
                s = -1

            if len(line) > 1:
                element[i][s] = line

    return element


def plot_data(data, outlier=False):
    c = sns.color_palette()
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
    plt.hlines(m, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m - 3*s, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m + 3*s, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    if z1[0] < -0.001:
        plt.plot(EP, p1(EP), color=c[2])
        print('EW slope: %.3f. Lower Teff' % z1[0])
    elif z1[0] > 0.001:
        plt.plot(EP, p1(EP), color=c[2])
        print('EW slope: %.3f. Higher Teff' % z1[0])
    else:
        plt.plot(EP, p1(EP), color=c[1])
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Excitation potential: $\chi$ [eV]')
    plt.ylabel('Abundance')

    plt.subplot(212)
    plt.plot(logRW, abund, '.w', label='Micro turbulence')
    plt.plot(logRW, abund, '.k')
    plt.hlines(m, min(logRW), max(logRW), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m - 3*s, min(logRW), max(logRW), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m + 3*s, min(logRW), max(logRW), colors=c[0], linestyles='--', linewidth=3)
    if z2[0] < -0.003:
        plt.plot(logRW, p2(logRW), color=c[2])
        print('RW slope: %.3f. Lower vt' % z2[0])
    if z2[0] > 0.003:
        plt.plot(logRW, p2(logRW), color=c[2])
        print('RW slope: %.3f. Higher vt' % z2[0])
    else:
        plt.plot(logRW, p2(logRW), color=c[1])
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Reduced EW: $\log RW$')
    plt.ylabel('Abundance')
    if outlier:
        idx = abs(data[:, 6]) > 3*s  # deviation larger than 3 sigma
        if True in idx:
            for i, w in enumerate(data[idx, 0]):
                ai = abund[idx][i]
                std = abs((ai-m)/s)
                print('Line at: %.3f is %.2f sigma away' % (w, std))
    return z1[0], z2[0]


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Force fit abundances with MOOG',
                                epilog='Happy spectroscopying :)')

    p.add_argument('teff', type=int, help='The effective temperature')
    p.add_argument('logg', type=float, help='The surface gravity')
    p.add_argument('feh', type=float, help='The metallicity')
    p.add_argument('vmicro', type=float, help='The micro turbulence')
    p.add_argument('-l', '--linelist', help='The linelist to be used')
    p.add_argument('-o', '--output', help='The output file with abundances', default='summary.out')
    p.add_argument('-p', '--plot', help='Enable plotting', default=True, type=bool)
    p.add_argument('-u', '--outlier', help='print 3 sigma outliers', default=False, action='store_true')

    args = p.parse_args()

    # Interpolate and transform models
    _create_model(args.teff, args.logg, args.feh, args.vmicro)

    # Update the batch file if a new linelist is used
    if args.linelist:
        _update_batch(args.linelist)

    # Run moog
    os.system('MOOGSILENT > /dev/null')

    # Prepare the data
    data = read_output(args.output)

    c1, c2 = plot_data(data[0], args.outlier)
    plt.tight_layout()
    plt.show()
