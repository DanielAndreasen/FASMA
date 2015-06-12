#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
import os
import numpy as np
from model_interpolation import save_model, interpolator
from pymoog import _get_model


def _run_moog(par='batch.par'):
    """Run MOOGSILENT with the given parameter file

    Raise an IOError if the parameter file does not exists
    """

    # if not os.path.exists(par):
        # raise IOError('The parameter file %s does not exists' % par)

    os.system('MOOGSILENT > zzz')
    os.system('rm -f zzz')

    # if par != 'batch.par':
    #     os.system('cp %s batch.par' % par)
    #     os.system('MOOGSILENT')
    #     os.system('rm -f batch.par')
    # else:
    #     os.system('MOOGSILENT %s' % par)


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


def fun_moog(x, par='batch.par', results='summary.out', fix_logg=False):
    """The 'function' that we should minimize

    :x: A tuple/list with values (teff, logg, [Fe/H], vt)
    :par: The parameter file (batch.par)
    :results: The summary file
    :returns: The slopes and abundances for the different elements
    """

    # Create an atmosphere model from input parameters
    teff, logg, feh, vt = x
    models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=feh)
    model = interpolator(models, teff=(teff, nt), logg=(logg, nl),
                         feh=(feh, nf))
    save_model(model, x)

    # Run MOOG and get the slopes and abundaces
    _run_moog(par=par)
    EPs, RWs, abundances = _read_moog(fname=results)
    if len(abundances) == 2:
        res = EPs[0]**2 + RWs[0]**2 + np.diff(abundances)[0]**2
        return res, EPs[0], RWs[0], abundances


def fun_moog_fortran(x, par='batch.par', results='summary.out', fix_logg=False):
    """The 'function' that we should minimize

    :x: A tuple/list with values (teff, logg, [Fe/H], vt)
    :par: The parameter file (batch.par)
    :results: The summary file
    :returns: The slopes and abundances for the different elements
    """

    import os
    # Create an atmosphere model from input parameters
    teff, logg, feh, vt = x
    p = '/home/daniel/Software/SPECPAR/interpol_models/'
    os.system('echo %i %s %s | %sintermod.e > z1' % (teff, logg, feh, p))
    os.system('echo %s | %stransform.e > z2' % (vt, p))
    os.system('rm -f z1 z2')

    # Run MOOG and get the slopes and abundaces
    _run_moog(par=par)
    EPs, RWs, abundances = _read_moog(fname=results)
    if len(abundances) == 2:
        res = EPs[0]**2 + RWs[0]**2 + np.diff(abundances)[0]**2
        return res, EPs[0], RWs[0], abundances
