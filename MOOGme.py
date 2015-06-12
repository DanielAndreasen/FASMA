#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import numpy as np
import argparse
import seaborn as sns
sns.set_style('dark')
sns.set_context('talk')
import matplotlib.pyplot as plt



def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Plot fits file for ARES. Be'
                                     ' careful with large files')


    parser.add_argument('linelist', help='Input line list')
    parser.add_argument('-p', '--parameter file',
            help='The parameter file (default: batch.par)',
            default='batch.par')
    parser.add_argument('-m', '--model',
            help='Model atmosphere',
            default='K95',
            choices=['K95', 'Kn', 'M', 'P'])
    parser.add_argument('-i', '--initial',
            help='Initial conditions (Teff, logg, [Fe/H], vt)',
            nargs='+',
            type=float,
            default=[5777, 4.44, 0.00, 1.00])
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
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = _parser()
    print(args)
