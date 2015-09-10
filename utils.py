#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import os
from model_interpolation import interpolator, save_model
import numpy as np
from glob import glob

K95 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                32000, 33000, 34000, 35000, 3500, 36000, 37000, 38000, 39000),
       'feh': (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
               0.1, 0.2, 0.3, 0.5, 1.0),
       'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}


def _get_model(teff, logg, feh, type='kurucz95'):
    """
    Find the names of the closest grid points for a given effective
    temperature, surface gravity, and iron abundance (proxy for metallicity).

    Inputs
    -----

    teff  :   The effective temperature(K) for the model atmosphere
    logg  :   The surface gravity (logarithmic in cgs) for the model atmosphere
    feh   :   The metallicity for the model atmosphere
    type  :   The type of atmosphere models to use. Currently only Kurucz from
              '95 is supported.

    output
    ------
    models      : List with path to 8 models the two closest in each parameter
                  space (2x2x2)
    teff_model  : The two closest effective temperatures in the grid
    logg_model  : The two closest surface gravities in the grid
    feh_model   : The two closest metallicities in the grid

    The last three return values are used for the interpolation to do some
    mapping. If only the paths to the models are needed, do not pay attention
    to them.

    """

    # Using the correct model atmosphere
    if type == 'kurucz95':
        grid = K95
    else:
        raise NotImplementedError('You request for type: %s is not available' % type)

    # Checking for bounds in Teff, logg, and [Fe/H]
    if (teff < grid['teff'][0]) or (teff > grid['teff'][-1]):
        raise ValueError('Teff out of bounds: %s' % teff)
    if (logg < grid['logg'][0]) or (logg > grid['logg'][-1]):
        raise ValueError('logg out of bounds: %s' % logg)
    if (feh < grid['feh'][0]) or (feh > grid['feh'][-1]):
        raise ValueError('[Fe/H] out of bounds: %s' % feh)

    if teff in grid['teff']:
        for i, Ti in enumerate(grid['teff']):
            if Ti - teff == 0:
                break
        teff_model = [grid['teff'][i], grid['teff'][i+1]]
    else:
        for i, Ti in enumerate(grid['teff']):
            if Ti - teff > 0:
                break
        teff_model = [grid['teff'][i-1], grid['teff'][i]]

    if logg in grid['logg']:
        for i, li in enumerate(grid['logg']):
            if li - logg == 0:
                break
        try:  # TODO: Do this everywhere!!!
            logg_model = [grid['logg'][i], grid['logg'][i+1]]
        except IndexError:
            logg_model = [grid['logg'][i-1], grid['logg'][i]]
    else:
        for i, li in enumerate(grid['logg']):
            if li - logg > 0:
                break
        logg_model = [grid['logg'][i-1], grid['logg'][i]]

    if feh in grid['feh']:
        for i, fi in enumerate(grid['feh']):
            if fi - feh == 0:
                break
        feh_model = [grid['feh'][i], grid['feh'][i+1]]
    else:
        for i, fi in enumerate(grid['feh']):
            if fi - feh > 0:
                break
        feh_model = [grid['feh'][i-1], grid['feh'][i]]

    if not os.path.isdir('models'):
        raise IOError('The models have to be inside a folder called models.')

    name = lambda t, g, s, f: 'models/kurucz95/%s%s/%ig%i.%s%s.gz' % (s, f, t, g*10, s, f)
    models = []
    for teff_m in teff_model:
        for logg_m in logg_model:
            for feh_m in feh_model:
                if feh_m >= 0:
                    feh_m = str(feh_m).replace('.', '')
                    fout = name(teff_m, logg_m, 'p', feh_m)
                else:
                    feh_m = str(feh_m).replace('.', '').replace('-', '')
                    fout = name(teff_m, logg_m, 'm', feh_m)
                models.append(fout)

    return models, teff_model, logg_model, feh_model


def _update_par(atmosphere_model='out.atm', line_list='linelist.moog', **kwargs):
    """Update the parameter file (batch.par) with new linelists, atmosphere
    models, or others.

    Inputs
    -----
    atmosphere_model    :   Location of your model atmosphere file
    line_list           :   Location of your line list

    Additional keyword arguments
    ----------------------------
    These additional keyword arguments allow the user to have full control
    over what is put into the MOOG input file. The default values are:

    terminal        'x11'
    atmosphere      1
    molecules       2
    trudamp         1
    lines           1
    flux/int        1
    damping         2
    units           0
    iraf            0
    plot            2
    obspectrum      1       Unless obspectrum is provided to the function.
    opacit          0
    freeform        0
    strong          0       Unless a strong lines list is provided.
    plotpars        1       0.75 Gaussian smoothing by default. Show full
                            synthesized spectral range with y:[0, 1.2]
    histogram       0
    synlimits               Defaults to the wavelength range provided and
                            the given wavelength step size, and the delta
                            defaults to the wavelength step size.

    Outputs
    -------
    And updated parameter file
    """

    # Path checks for input files
    if not os.path.exists(line_list):
        raise IOError('Line list file "%s" could not be found.' % (line_list))

    default_kwargs = {
        'atmosphere': 1,
        'molecules':  1,
        'trudamp':    1,  # Sure, why not? It's a black art anyway!
        'lines':      1,
        'terminal':   'x11',
        'flux/int':   0,
        'damping':    2,
        'units':      0,
        'iraf':       0,
        'plot':       0,
        'obspectrum': 0,
        'opacit':     0,
        'freeform':   0,
        'strong':     0,
        'summary':    'summary.out'
        }

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.iteritems():
        if key not in kwargs.keys():
            kwargs[key] = value
    # Generate a MOOG-compatible run file
    moog_filename = 'batch.par'

    moog_contents = "abfind\n"\
                    "terminal       %s\n"\
                    "model_in       '%s'\n"\
                    "summary_out    '%s'\n"\
                    "standard_out   '%s'\n"\
                    "lines_in       '%s'\n" % (kwargs['terminal'], atmosphere_model,
                                               kwargs['summary'], 'result.out', line_list)

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,plot,opacit,freeform,obspectrum,histogram,'\
               'synlimits'.split(',')
    if 'plotpars' in kwargs:
        if kwargs['plotpars'] != 0:
            settings.append('plotpars')

    for setting in settings:
        if setting in kwargs:
            moog_contents += "%s %s\n" % (setting + ' ' * (14 - len(setting)), kwargs[setting])

    with open(moog_filename, 'w') as moog:
        moog.writelines(moog_contents)


def _run_moog(par='batch.par'):
    """Run MOOGSILENT with the given parameter file
    """
    os.system('MOOGSILENT > /dev/null')


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


def fun_moog(x, par='batch.par', results='summary.out'):
    """The 'function' that we should minimize

    :x: A tuple/list with values (teff, logg, [Fe/H], vt)
    :par: The parameter file (batch.par)
    :results: The summary file
    :returns: The slopes and abundances for the different elements
    """

    # Create an atmosphere model from input parameters
    teff, logg, feh, _ = x
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


def fun_moog_fortran(x, par='batch.par', results='summary.out'):
    """The 'function' that we should minimize

    :x: A tuple/list with values (teff, logg, [Fe/H], vt)
    :par: The parameter file (batch.par)
    :results: The summary file
    :returns: The slopes and abundances for the different elements
    """

    # Create an atmosphere model from input parameters
    teff, logg, feh, vt = x
    p = '/home/daniel/Software/SPECPAR/interpol_models/'
    os.system('echo %i %s %s | %sintermod.e > /dev/null' % (teff, logg, feh, p))
    os.system('echo %s | %stransform.e > /dev/null' % (vt, p))

    # Run MOOG and get the slopes and abundaces
    _run_moog(par=par)
    EPs, RWs, abundances = _read_moog(fname=results)
    if len(abundances) == 2:
        res = EPs[0]**2 + RWs[0]**2 + np.diff(abundances)[0]**2
        return res, EPs[0], RWs[0], abundances


def readmoog(output):
    """Read the output file from MOOG"""

    nelements = 1
    readdata = False
    Fe1Lines = []
    Fe2Lines = []
    with open(output, 'r') as lines:
        for line in lines:
            if 'Teff' in line:  # Get the atmospheric parameters
                line = line.split()
                teff = int(line[1])
                logg = float(line[4])
                vt = float(line[6])
                feh = float(line[-1].split('=')[-1])
            elif '#lines' in line and nelements == 1:  # Statistics on FeI
                readdata = False
                line = line.split()
                nfe1 = int(line[-1])
                fe1 = float(line[3])
                sigfe1 = float(line[7])
            elif '#lines' in line and nelements == 2:  # Statistics on FeII
                readdata = False
                line = line.split()
                nfe2 = int(line[-1])
                fe2 = float(line[3])
                sigfe2 = float(line[7])
            elif 'E.P.' in line and nelements == 1:  # We only want information from FeI
                line = line.split()
                slopeEP = float(line[4])
            elif 'R.W.' in line and nelements == 1:  # We only want information from FeI
                line = line.split()
                nelements += 1  # Done with this element, move to next one
                try:
                    slopeRW = float(line[4])
                except ValueError:
                    slopeRW = False
            else:
                if line.startswith('wavelength'):
                    readdata = True
                    continue
            if readdata:
                content = map(float, filter(None, line.split(' ')))
                if nelements == 1:
                    Fe1Lines.append(content)
                else:
                    Fe2Lines.append(content)

    # Store the line information in numpy arrays because lists are not for science!
    linesFe1 = np.zeros((len(Fe1Lines), 7))
    linesFe2 = np.zeros((len(Fe2Lines), 7))
    for i, f1 in enumerate(Fe1Lines):
        linesFe1[i, 0] = f1[0]
        linesFe1[i, 1] = f1[1]
        linesFe1[i, 2] = f1[2]
        linesFe1[i, 3] = f1[3]
        linesFe1[i, 4] = f1[4]
        linesFe1[i, 5] = f1[5]
        linesFe1[i, 6] = f1[6]
    for i, f2 in enumerate(Fe2Lines):
        linesFe2[i, 0] = f2[0]
        linesFe2[i, 1] = f2[1]
        linesFe2[i, 2] = f2[2]
        linesFe2[i, 3] = f2[3]
        linesFe2[i, 4] = f2[4]
        linesFe2[i, 5] = f2[5]
        linesFe2[i, 6] = f2[6]

    # If We don't have any RW slope, calculate it manually
    if not slopeRW:
        slopeRW, _ = np.polyfit(linesFe1[:, 4], linesFe1[:, 5], 1)
    sigfe1 = sigfe1 / np.sqrt(nfe1)
    sigfe2 = sigfe2 / np.sqrt(nfe2)
    return teff, logg, vt, feh, fe1-7.47, sigfe1, fe2-7.47, sigfe2, slopeEP, slopeRW, linesFe1, linesFe2


def _slopeSigma(x, y):
    """Sigma on a slope after fitting a straight line"""
    N = len(x)
    sxx = np.sum((x-np.mean(x))**2)
    a, b = np.polyfit(x, y, 1)
    chi2 = np.sum((y - a*x-b)**2)
    return np.sqrt(chi2/((N-2)*sxx))


def error(linelist):
    """linelist to give error estimation on"""
    # Find the output file and read the current state of it
    if os.path.isfile('results/%s.out' % linelist):
        summary = readmoog('results/%s.out' % linelist)
    else:
        summary = readmoog('results/%s.NC.out' % linelist)
    # Read the correct output file (error_summary.out).
    _update_par(line_list='linelist/%s' % linelist, summary='error_summary.out')

    # Prepare the different things we need
    teff, logg, vt, feh = summary[0:4]
    Fe1, Fe2 = summary[-2], summary[-1]
    sigmafe1 = summary[5]
    sigmafe2 = summary[7]

    siga1 = _slopeSigma(Fe1[:, 4], Fe1[:, 5])
    siga2 = _slopeSigma(Fe1[:, 1], Fe1[:, 5])

    # Error om microturbulence
    fun_moog_fortran((teff, logg, feh, vt+0.1), results='error_summary.out')
    sumvt = readmoog('error_summary.out')
    slopeEP, slopeRW = sumvt[8], sumvt[9]
    if slopeRW == 0:
        errormicro = abs(siga1/0.001) * 0.10
    else:
        errormicro = abs(siga1/slopeRW) * 0.10

    # Contribution to [Fe/H]
    deltafe1micro = abs((errormicro/0.10) * (sumvt[4]-feh))

    # Error on Teff
    slopes = errormicro/0.10 * slopeEP
    errorslopeEP = np.hypot(slopes, siga2)
    fun_moog_fortran((teff+100, logg, feh, vt), results='error_summary.out')
    sumteff = readmoog('error_summary.out')

    errorteff = abs(errorslopeEP/sumteff[8]) * 100
    # Contribution to [Fe/H]
    deltafe1teff = abs((errorteff/100) * (sumteff[4]-feh))

    # Error on logg
    fe2error = abs(errorteff/100 * sumteff[6]-feh)
    sigmafe2total = np.hypot(sigmafe2, fe2error)
    fun_moog_fortran((teff, logg-0.20, feh, vt), results='error_summary.out')
    sumlogg = readmoog('error_summary.out')
    errorlogg = abs(sigmafe2total/(sumlogg[6]-feh)*0.20)

    # Error on [Fe/H]
    errorfeh = np.sqrt(sigmafe1**2 + deltafe1teff**2 + deltafe1micro**2)

    errorteff = int(errorteff)
    errorlogg = round(errorlogg, 2)
    errorfeh = round(errorfeh, 2)
    errormicro = round(errormicro, 2)

    os.remove('error_summary.out')
    return teff, errorteff, logg, errorlogg, feh, errorfeh, vt, errormicro