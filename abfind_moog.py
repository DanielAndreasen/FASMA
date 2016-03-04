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
from utils import Readmoog

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


def _update_batch(linelist=False):
    """Update the batch.par file for MOOG"""
    with open('batch.par', 'r') as lines:
        lines_tmp = []
        for line in lines:
            if ('lines_in' in line) and linelist:
                line = line.split('      ')
                line[1] = "'" + linelist + "'\n"
                line = '       '.join(line)
            if 'summary_out' in line:
                line = line.split()
                line[1] = "'%s'\n" % 'summary.out'
                line = '    '.join(line)
            lines_tmp.append(line)

    with open('batch.par', 'w') as lines:
        lines.write(''.join(lines_tmp))


def plot_data(data, outlier=False, version=2014):
    idx = 1 if version > 2013 else 0
    c = sns.color_palette()
    EP = data[:, 1+idx]
    logRW = data[:, 4+idx]
    abund = data[:, 5+idx]
    m = np.mean(abund)
    s = 3*np.std(abund)

    z1 = np.polyfit(EP, abund, 1)
    p1 = np.poly1d(z1)

    z2 = np.polyfit(logRW, abund, 1)
    p2 = np.poly1d(z2)

    plt.figure(figsize=(8, 9))
    plt.subplot(211)
    plt.plot(EP, abund, '.w', label='Temperature')
    plt.plot(EP, abund, '.k')
    plt.hlines(m, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m - s, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m + s, min(EP), max(EP), colors=c[0], linestyles='--', linewidth=3)
    if z1[0] < -0.001:
        plt.plot(EP, p1(EP), color=c[2], lw=3)
        print('EW slope: %.3f. Lower Teff' % z1[0])
    elif z1[0] > 0.001:
        plt.plot(EP, p1(EP), color=c[2], lw=3)
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
    plt.hlines(m - s, min(logRW), max(logRW), colors=c[0], linestyles='--', linewidth=3)
    plt.hlines(m + s, min(logRW), max(logRW), colors=c[0], linestyles='--', linewidth=3)
    if z2[0] < -0.003:
        plt.plot(logRW, p2(logRW), color=c[2], lw=3)
        print('RW slope: %.3f. Lower vt' % z2[0])
    elif z2[0] > 0.003:
        plt.plot(logRW, p2(logRW), color=c[2], lw=3)
        print('RW slope: %.3f. Higher vt' % z2[0])
    else:
        plt.plot(logRW, p2(logRW), color=c[1])
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Reduced EW: $\log RW$')
    plt.ylabel('Abundance')
    if outlier:
        indc = abs(data[:, 6+idx]) > s  # deviation larger than 3 sigma
        if True in indc:
            for i, w in enumerate(data[indc, 0]):
                ai = abund[indc][i]
                std = abs((ai-m)/s)
                print('Line at: %.3f is %.2f sigma away' % (w, std))
    return z1[0], z2[0]


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Force fit abundances with MOOG', epilog='Happy spectroscopying :)')
    p.add_argument('teff', type=int, help='The effective temperature')
    p.add_argument('logg', type=float, help='The surface gravity')
    p.add_argument('feh', type=float, help='The metallicity')
    p.add_argument('vmicro', type=float, help='The micro turbulence')
    p.add_argument('-l', '--linelist', help='The linelist to be used', default=False)
    p.add_argument('-p', '--plot', help='Enable plotting', default=True, action='store_false')
    p.add_argument('-u', '--outlier', help='print 3 sigma outliers', default=False, action='store_true')
    p.add_argument('-v', '--version', help='MOOG version', choices=['2013', '2014'], default='2014')
    args = p.parse_args()

    # Interpolate and transform models
    _create_model(args.teff, args.logg, args.feh, args.vmicro)

    # Update the batch file if a new linelist is used
    _update_batch(args.linelist)

    # Run moog
    os.system('MOOGSILENT > /dev/null')

    # Prepare the data
    m = Readmoog(fname='summary.out', version=int(args.version))
    _, _, _, _, _, _, data, _ = m.fe_statistics()
    c1, c2 = plot_data(data, args.outlier, int(args.version))
    if args.plot:
        plt.tight_layout()
        plt.show()
