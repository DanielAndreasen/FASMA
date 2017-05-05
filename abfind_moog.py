#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from utils import Readmoog
from interpolation import interpolator
import statsmodels.formula.api as sm
import os
import argparse


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
    wave = data[:, 0]
    EP = data[:, 1+idx]
    logRW = data[:, 4+idx]
    abund = data[:, 5+idx]
    sig = data[:, 6+idx]
    m = np.mean(abund)
    s = 3*np.std(abund)

    data = {'EP': EP, 'logRW': logRW, 'abund': abund, 'sig': sig, 'wave': wave}
    w = np.ones(len(data['abund']))

    z1 = sm.wls('abund ~ EP', data=data, weights=w).fit().params
    z1 = [z1[1], z1[0]]
    p1 = np.poly1d(z1)

    z2 = sm.wls('abund ~ logRW', data=data, weights=w).fit().params
    z2 = [z2[1], z2[0]]
    p2 = np.poly1d(z2)

    plt.figure(figsize=(8, 9))
    plt.subplot(211)
    plt.plot(EP, abund, '.w', label='Temperature')
    plt.plot(EP, abund, '.k')
    plt.hlines(m, min(EP), max(EP), colors='C0', linestyles='--', linewidth=3)
    plt.hlines(m - s, min(EP), max(EP), colors='C0', linestyles='--', linewidth=3)
    plt.hlines(m + s, min(EP), max(EP), colors='C0', linestyles='--', linewidth=3)
    if z1[0] < -0.001:
        plt.plot(EP, p1(EP), color='C3', lw=3)
        print('EW slope: %.3f. Lower Teff' % z1[0])
    elif z1[0] > 0.001:
        plt.plot(EP, p1(EP), color='C3', lw=3)
        print('EW slope: %.3f. Higher Teff' % z1[0])
    else:
        plt.plot(EP, p1(EP), color='C2')
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Excitation potential: $\chi$ [eV]')
    plt.ylabel('Abundance')
    y1, y2 = plt.ylim()
    plt.ylim(0.9*y1, 1.1*y2)

    plt.subplot(212)
    plt.plot(logRW, abund, '.w', label='Micro turbulence')
    plt.plot(logRW, abund, '.k')
    plt.hlines(m, min(logRW), max(logRW), colors='C0', linestyles='--', linewidth=3)
    plt.hlines(m - s, min(logRW), max(logRW), colors='C0', linestyles='--', linewidth=3)
    plt.hlines(m + s, min(logRW), max(logRW), colors='C0', linestyles='--', linewidth=3)
    if z2[0] < -0.003:
        plt.plot(logRW, p2(logRW), color='C3', lw=3)
        print('RW slope: %.3f. Lower vt' % z2[0])
    elif z2[0] > 0.003:
        plt.plot(logRW, p2(logRW), color='C3', lw=3)
        print('RW slope: %.3f. Higher vt' % z2[0])
    else:
        plt.plot(logRW, p2(logRW), color='C2')
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'Reduced EW: $\log RW$')
    plt.ylabel('Abundance')
    plt.ylim(0.9*y1, 1.1*y2)
    if outlier:
        indc = abs(data['sig']) > s  # deviation larger than 3 sigma
        if True in indc:
            for ai, w in zip(data['abund'][indc], data['wave'][indc]):
                std = 3*abs((ai-m)/s)
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
    x = (args.teff, args.logg, args.feh, args.vmicro)
    interpolator(params=x)

    # Update the batch file if a new linelist is used
    _update_batch(args.linelist)

    # Run moog
    os.system('MOOGSILENT > /dev/null')

    # Prepare the data
    m = Readmoog(params=x, fname='summary.out', version=int(args.version))
    f1, _, f2, _, _, _, fe1, fe2 = m.fe_statistics()
    c1, c2 = plot_data(fe1, args.outlier, int(args.version))
    c1, c2 = plot_data(fe2, args.outlier, int(args.version))
    f12 = f1-f2
    if f12 < 0.0:
        print('Fe diff: %.2f. Lower logg' % f12)
    elif f12 > 0.0:
        print('Fe diff: %.2f. Higher logg' % f12)
    if args.plot:
        plt.tight_layout()
        plt.show()
