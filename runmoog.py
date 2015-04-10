#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import os


def _run_moog(par='batch.par'):
    """Run MOOGSILENT with the given parameter file

    Raise an IOError if the parameter file does not exists
    """

    if not os.path.exists(par):
        raise IOError('The parameter file %s does not exists' % par)

    if par != 'batch.par':
        os.system('cp %s batch.par' % par)
        os.system('MOOGSILENT %s' % par)
        os.system('rm -f batch.par')
    else:
        os.system('MOOGSILENT %s' % par)


def _read_moog(fname='summary.out'):
    """Read the slopes from the summary.out and return them

    :fname: From the summary_out
    :returns: A tuple of the slopes and the average abundances for
    different elements
    """
    EP_slopes = []
    RW_slopes = []
    abundances = []
    with open(fname, 'r') as lines:
        for line in lines:
            # Get the EP slope
            if line.startswith('E.P'):
                line = filter(None, line.split('slope =')[1].split(' '))
                EP_slopes.append(float(line[0]))
            # Get the reduced EW slope
            elif line.startswith('R.W'):
                line = filter(None, line.split('slope =')[1].split(' '))
                RW_slopes.append(float(line[0]))
            # Get the average abundance
            elif line.startswith('average abundance'):
                line = filter(None, line.split('abundance =')[1].split(' '))
                abundances.append(float(line[0]))

    return EP_slopes, RW_slopes, abundances


def fun_moog(par='batch.par', results='summary.out'):
    """The 'function' that we should minimize

    :par: The parameter file (batch.par)
    :results: The summary file
    :returns: The slopes and abundances for the different elements

    """
    _run_moog(par=par)
    return _read_moog(fname=results)

if __name__ == '__main__':
    print fun_moog()
