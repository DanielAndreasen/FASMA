#!/usr/bin/env python
# -*- coding: utf8 -*-

from __future__ import division
import os
from interpolation import save_model, interpolator
import numpy as np
from glob import glob

kurucz95 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                32000, 33000, 34000, 35000, 3500, 36000, 37000, 38000, 39000),
       'feh': (-3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
               0.1, 0.2, 0.3, 0.5, 1.0),
       'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}

kurucz08 = {'teff': (3750, 4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000,
                6250, 6500, 6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500,
                8750, 9000, 9250, 9500, 9750, 10000, 10250, 10500, 10750,
                11000, 11250, 11500, 11750, 12000, 12250, 12500, 12750, 13000,
                14000, 15000, 16000, 17000, 18000, 19000, 20000, 21000, 22000,
                23000, 24000, 25000, 26000, 27000, 28000, 29000, 30000, 31000,
                32000, 33000, 34000, 35000, 3500, 36000, 37000, 38000, 39000),
       'feh': (-4.0, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, -0.3, -0.2, -0.1, 0.0,
               0.1, 0.2, 0.3, 0.5, 1.0),
       'logg': (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0)}


def _get_model(teff, logg, feh, atmtype='kurucz95'):
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

    def _model_path(atmpath, feh, teff, logg):
        name = 'models/%s/' % atmpath
        if feh < 0:
            name += 'm%s/' % str(abs(feh)).replace('.', '')
        else:
            name += 'p%s/' % str(abs(feh)).replace('.', '')
        name += '%ig%s.' % (teff, str(logg).replace('.', ''))
        if feh < 0:
            name += 'm%s.gz' % str(abs(feh)).replace('.', '')
        else:
            name += 'p%s.gz' % str(abs(feh)).replace('.', '')
        return name

    # Using the correct model atmosphere
    atmmodels = {'kurucz95': [kurucz95, 'kurucz95'], 'kurucz08': [kurucz08, 'kurucz08']}
    if atmtype in atmmodels.keys():
        grid = atmmodels[atmtype][0]
    else:
        raise NotImplementedError('You request for atmospheric models: %s is not available' % atmtype)

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
        try:
            # teff_model = [grid['teff'][i], grid['teff'][i+1]]
            teff_model = [grid['teff'][i-1], grid['teff'][i], grid['teff'][i+1], grid['teff'][i+2]]
        except IndexError:
            teff_model = [grid['teff'][i-1], grid['teff'][i]]
    else:
        for i, Ti in enumerate(grid['teff']):
            if Ti - teff > 0:
                break
        try:
            # teff_model = [grid['teff'][i-1], grid['teff'][i]]
            teff_model = [grid['teff'][i-1], grid['teff'][i], grid['teff'][i+1], grid['teff'][i+2]]
        except IndexError:
            teff_model = [grid['teff'][i], grid['teff'][i+1]]

    if logg in grid['logg']:
        for i, li in enumerate(grid['logg']):
            if li - logg == 0:
                break
        try:
            logg_model = [grid['logg'][i], grid['logg'][i+1]]
        except IndexError:
            logg_model = [grid['logg'][i-1], grid['logg'][i]]
    else:
        for i, li in enumerate(grid['logg']):
            if li - logg > 0:
                break
        try:
            logg_model = [grid['logg'][i-1], grid['logg'][i]]
        except IndexError:
            logg_model = [grid['logg'][i], grid['logg'][i+1]]

    if feh in grid['feh']:
        for i, fi in enumerate(grid['feh']):
            if fi - feh == 0:
                break
        try:
            feh_model = [grid['feh'][i], grid['feh'][i+1]]
        except IndexError:
            feh_model = [grid['feh'][i-1], grid['feh'][i]]
    else:
        for i, fi in enumerate(grid['feh']):
            if fi - feh > 0:
                break
        try:
            feh_model = [grid['feh'][i-1], grid['feh'][i]]
        except IndexError:
            feh_model = [grid['feh'][i], grid['feh'][i+1]]

    if not os.path.isdir('models'):
        raise IOError('The models have to be inside a folder called models.')

    models = []
    for teff_m in teff_model:
        for logg_m in logg_model:
            for feh_m in feh_model:
                models.append(_model_path(atmmodels[atmtype][1], feh_m, teff_m, logg_m))

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

    with open('batch.par', 'w') as moog:
        moog.writelines(moog_contents)


def _update_par_synth(start_wave, end_wave, line_list='linelist.moog', atmosphere_model='out.atm', **kwargs):
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
                            synthesized spectral range with y:[0, 1.2]vsini
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
        'molecules':  2,
        'lines':      1,
        'terminal':   'x11',
        'flux/int':   0,
        'damping':    2,
        'obspectrum': 0,
        'model_in':     'out.atm',
        'smoothed_out': 'smooth.out',
        'summary':      'summary.out'
        }

    # Fill the keyword arguments with the defaults if they don't exist already
    for key, value in default_kwargs.iteritems():
        if key not in kwargs.keys():
            kwargs[key] = value
    # Generate a MOOG-compatible run file

    moog_contents = "synth\n"\
                    "terminal          %s\n"\
                    "model_in          '%s'\n"\
                    "observed_in       '%s'\n"\
                    "summary_out       '%s'\n"\
                    "smoothed_out      'smooth.out'\n"\
                    "standard_out      'result.out'\n"\
                    "lines_in          '%s'\n"\
                    "plot              1\n"\
                    "synlimits\n"\
                    "      %s      %s       %s      %s\n"\
                    "plotpars          %s\n"\
                    "      %s      %s       0.5      1.05\n"\
                    "      0.0     0.0      0.0       0.0\n"\
                    "      g       %s       %s       %s       %s       %s\n" % (kwargs['terminal'], atmosphere_model,  kwargs['obfile'], kwargs['summary'],
                                                                               line_list, start_wave, end_wave, kwargs['step_wave'], kwargs['step_flux'],
                                                                               kwargs['plotpars'], start_wave, end_wave, kwargs['resolution'], kwargs['vsini'],
                                                                               kwargs['limb'], kwargs['vmac'], kwargs['lorentz'])

    settings = 'atmosphere,molecules,trudamp,lines,strong,flux/int,damping,'\
               'units,iraf,opacity,freeform,obspectrum,histogram,'\
               'synlimits'.split(',')

    for setting in settings:
        if setting in kwargs:
            moog_contents += "%s      %s\n" % (setting + ' ' * (14 - len(setting)), kwargs[setting])

    with open('batch.par', 'w') as moog:
        moog.writelines(moog_contents)


def _run_moog(par='batch.par', driver='abfind'):
    """Run MOOGSILENT with the given parameter file"""
    if driver == 'abfind':
        os.system('MOOGSILENT > /dev/null')
        # os.system('/home/daniel/Software/moogjul2014/MOOGSILENT > /dev/null')
    elif driver == 'synth':
        with open('stupid.tmp', 'w') as f:
            f.writelines('batch.par\nq')
        os.system('MOOGSILENT < stupid.tmp > /dev/null')
        os.remove('stupid.tmp')


def _read_smooth(fname='smooth.out'):
    """Read the synthetic spectrum from the summary.out and return them
    :fname: From the summary_out
    :returns: wavelength and flux (in that order)
    """
    wavelength, flux = np.loadtxt(fname, skiprows=2, usecols=(0, 1), unpack=True)
    return wavelength, flux


def fun_moog(x, par='batch.par', results='summary.out', weights='null',
             driver='abfind', version=2013):
    """Run MOOG and return slopes for abfind mode.

    :x: A tuple/list with values (teff, logg, [Fe/H], vt)
    :par: The parameter file (batch.par)
    :results: The summary file
    :weights: The weights to be used in the slope calculation
    :driver: Which driver to use when running MOOG
    :version: The version of MOOG
    :returns: The slopes and abundances for the different elements
    """

    # Create an atmosphere model from input parameters
    teff, logg, feh, _ = x
    models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=feh)
    model = interpolator(models, teff=(teff, nt), logg=(logg, nl), feh=(feh, nf))
    save_model(model, x)

    # Run MOOG and get the slopes and abundaces
    _run_moog(par=par, driver=driver)
    if driver == 'abfind':
        m = Readmoog(fname=results, version=version)
        _, _, _, _, _, _, data, _ = m.fe_statistics()
        if version > 2013:
            EPs, _ = slope((data[:,2], data[:,6]), weights=weights)
            RWs, _ = slope((data[:,5], data[:,6]), weights=weights)
        else:
            EPs, _ = slope((data[:,1], data[:,5]), weights=weights)
            RWs, _ = slope((data[:,4], data[:,5]), weights=weights)
        m = Readmoog(fname=results, version=version)
        fe1, _, fe2, _, _, _, _, _ = m.fe_statistics()
        abundances = [fe1+7.47, fe2+7.47]
        res = EPs**2 + RWs**2 + np.diff(abundances)[0]**2
        return res, EPs, RWs, abundances


class Readmoog:

    def __init__(self, fname='summary.out', version=2013):
        self.fname = fname
        self.nelements = 1
        self.idx = 1 if version > 2013 else 0
        self.version = version
        with open(self.fname, 'r') as f:
            self.lines = f.readlines()

    def parameters(self):
        """Get the atmospheric parameters"""
        for line in self.lines:
            if 'Teff' in line:
                break
        line = line.split()
        self.teff = int(line[1])
        self.logg = float(line[4])
        self.vt = float(line[6])
        self.feh = float(line[-1].split('=')[-1])
        self.params = self.teff, self.logg, self.feh, self.vt
        return self.params


    def fe_statistics(self):
        """Get statistics on Fe lines"""
        self.readdata = False
        self.Fe1Lines = []
        self.Fe2Lines = []
        for line in self.lines:
            if '#lines' in line and self.nelements == 1:  # Statistics on FeI
                line = line.split()
                self.readdata = False
                self.nfe1 = int(line[-1])
                self.fe1 = float(line[3])
                self.sigfe1 = float(line[7])
            elif '#lines' in line and self.nelements == 2:  # Statistics on FeII
                line = line.split()
                self.readdata = False
                self.nfe2 = int(line[-1])
                self.fe2 = float(line[3])
                self.sigfe2 = float(line[7])
            elif 'E.P.' in line and self.nelements == 1:  # We only want information from FeI
                line = line.split()
                try:
                    self.slopeEP = float(line[4])
                except ValueError:
                    self.slopeEP = False
            elif 'R.W.' in line and self.nelements == 1:  # We only want information from FeI
                line = line.split()
                self.nelements += 1  # Done with this element, move to next one
                try:
                    self.slopeRW = float(line[4])
                except ValueError:
                    self.slopeRW = False
            else:
                if line.startswith('wavelength'):
                    self.readdata = True
                    continue
            if self.readdata:
                content = map(float, filter(None, line.split(' ')))
                if self.nelements == 1:
                    self.Fe1Lines.append(content)
                else:
                    self.Fe2Lines.append(content)

        # Store the line information in numpy arrays because lists are not for science!
        self.linesFe1 = np.zeros((len(self.Fe1Lines), 7+self.idx))
        self.linesFe2 = np.zeros((len(self.Fe2Lines), 7+self.idx))
        for i, f1 in enumerate(self.Fe1Lines):
            self.linesFe1[i, 0] = f1[0]
            self.linesFe1[i, 1] = f1[1]
            self.linesFe1[i, 2] = f1[2]
            self.linesFe1[i, 3] = f1[3]
            self.linesFe1[i, 4] = f1[4]
            self.linesFe1[i, 5] = f1[5]
            self.linesFe1[i, 6] = f1[6]
            if self.version > 2013:
                self.linesFe1[i, 7] = f1[7]
        for i, f2 in enumerate(self.Fe2Lines):
            self.linesFe2[i, 0] = f2[0]
            self.linesFe2[i, 1] = f2[1]
            self.linesFe2[i, 2] = f2[2]
            self.linesFe2[i, 3] = f2[3]
            self.linesFe2[i, 4] = f2[4]
            self.linesFe2[i, 5] = f2[5]
            self.linesFe2[i, 6] = f2[6]
            if self.version > 2013:
                self.linesFe2[i, 7] = f2[7]

        # If We don't have any RW slope, calculate it manually
        if not self.slopeRW:
            self.slopeRW, _ = np.polyfit(self.linesFe1[:, 4+self.idx], self.linesFe1[:, 5+self.idx], 1)
        if not self.slopeEP:
            self.slopeEP, _ = np.polyfit(self.linesFe1[:, 1+self.idx], self.linesFe1[:, 5+self.idx], 1)
        self.sigfe1 = self.sigfe1 / np.sqrt(self.nfe1)
        self.sigfe2 = self.sigfe2 / np.sqrt(self.nfe2)
        return self.fe1-7.47, self.sigfe1, self.fe2-7.47, self.sigfe2, self.slopeEP, self.slopeRW, self.linesFe1, self.linesFe2


    def elements(self):
        """Get the elements and abundances from the output file"""
        abundances = []
        element = []
        for line in self.lines:
            # Get the average abundance
            if line.startswith('average abundance'):
                line = filter(None, line.split('abundance =')[1].split(' '))
                abundances.append(float(line[0]))
              # Get element
            elif line.startswith('Abundance'):
                line = filter(None, line.split(' '))
                element.append(str(line[4])+str(line[5]))
        return element, abundances


def _slopeSigma(x, y, weights):
    """Sigma on a slope after fitting a straight line"""
    N = len(x)
    sxx = np.sum((x-np.mean(x))**2)
    a, b = np.polyfit(x, y, 1, w=weights)
    chi2 = np.sum((y - a*x-b)**2)
    return np.sqrt(chi2/((N-2)*sxx))


def error(linelist, converged, version=2013, weights='null'):
    """linelist to give error estimation on"""
    # Find the output file and read the current state of it
    idx = 1 if version > 2013 else 0
    if converged:
        m = Readmoog(fname='results/%s.out' % linelist, version=version)
        summary = m.fe_statistics()
    else:
        m = Readmoog(fname='results/%s.NC.out' % linelist, version=version)
        summary = m.fe_statistics()
    # Read the correct output file (error_summary.out).
    _update_par(line_list='linelist/%s' % linelist, summary='error_summary.out')
    data = summary[6]
    _, weights = slope((data[:,1+idx], data[:,5+idx]), weights=weights)

    # Prepare the different things we need
    teff, logg, feh, vt = m.parameters()
    Fe1, Fe2 = summary[-2], summary[-1]
    sigmafe1 = summary[1]
    sigmafe2 = summary[3]

    siga1 = _slopeSigma(Fe1[:, 4+idx], Fe1[:, 5+idx], weights=weights)
    siga2 = _slopeSigma(Fe1[:, 1+idx], Fe1[:, 5+idx], weights=weights)

    # Error om microturbulence
    fun_moog((teff, logg, feh, vt+0.1), results='error_summary.out', version=version)
    sumvt = Readmoog(fname='error_summary.out', version=version).fe_statistics()
    slopeEP, slopeRW = sumvt[4], sumvt[5]
    if slopeRW == 0:
        errormicro = abs(siga1/0.001) * 0.10
    else:
        errormicro = abs(siga1/slopeRW) * 0.10

    # Contribution to [Fe/H]
    deltafe1micro = abs((errormicro/0.10) * (sumvt[0]-feh))

    # Error on Teff
    slopes = errormicro/0.10 * slopeEP
    errorslopeEP = np.hypot(slopes, siga2)
    fun_moog((teff+100, logg, feh, vt), results='error_summary.out', version=version)
    sumteff = Readmoog(fname='error_summary.out', version=version).fe_statistics()

    errorteff = abs(errorslopeEP/sumteff[4]) * 100
    # Contribution to [Fe/H]
    deltafe1teff = abs((errorteff/100) * (sumteff[0]-feh))
    # Error on logg
    fe2error = abs(errorteff/100 * (sumteff[2]-feh))
    sigmafe2total = np.hypot(sigmafe2, fe2error)
    fun_moog((teff, logg-0.20, feh, vt), results='error_summary.out', version=version)
    sumlogg = Readmoog(fname='error_summary.out', version=version).fe_statistics()
    errorlogg = abs(sigmafe2total/(sumlogg[2]-feh)*0.20)

    # Error on [Fe/H]
    errorfeh = np.sqrt(sigmafe1**2 + deltafe1teff**2 + deltafe1micro**2)

    errorteff = int(errorteff)
    errorlogg = round(errorlogg, 2)
    errorfeh = round(errorfeh, 2)
    errormicro = round(errormicro, 2)

    os.remove('error_summary.out')
    return teff, errorteff, logg, errorlogg, feh, errorfeh, vt, errormicro


def slope(data, weights='null'):
    """Calculate the slope of a data set with weights"""
    import statsmodels.formula.api as sm
    weights = weights.lower()
    options = ['null', 'sigma', 'mad']
    if weights not in options:
        weights = None

    data = {'x': data[0], 'y': data[1]}
    fit = np.polyfit(data['x'], data['y'], 1)
    Y = np.poly1d(fit)(data['x'])
    dif = data['y'] - Y

    if not weights:
        w = 1/abs(dif)
        idx = np.isinf(w)
        w[~idx] /= w[~idx].max()
        w[idx] = 1
    if weights == 'null':
        w = np.ones(len(data['x']))
    elif weights == 'sigma':
        sig = np.std(dif)
        w = np.zeros(len(data['y'])) + 0.01
        mask3 = abs(data['y']-Y) < 3*sig
        w[mask3] = 0.10
        mask2 = abs(data['y'] - Y) < 2*sig
        w[mask2] = 0.25
        mask1 = abs(data['y'] - Y) < sig
        w[mask1] = 1.0
    elif weights == 'mad':
        mad = np.mean(np.absolute(dif - np.mean(dif, None)), None)
        w = np.zeros(len(data['y'])) + 0.01
        mask3 = abs(data['y']-Y) < 3*mad
        w[mask3] = 0.10
        mask2 = abs(data['y'] - Y) < 2*mad
        w[mask2] = 0.25
        mask1 = abs(data['y'] - Y) < mad
        w[mask1] = 1.0

    wls = sm.wls('y ~ x', data=data, weights=w).fit()
    return wls.params[1], w


def read_observations(wavelength, flux, start_synth, end_synth):
    """Read observed spectrum and return wavelength and flux"""
    wavelength_obs = wavelength[(wavelength >= start_synth) & (wavelength <= end_synth)]
    flux_obs = flux[(wavelength >= start_synth) & (wavelength <= end_synth)]
    return wavelength_obs, flux_obs


def interpol_synthetic(wavelength, flux, start_synth, end_synth):
    """Interpolation of the synthetic flux to the observed wavelength"""
    from scipy.interpolate import interp1d
    # The synthetic spectrum should be always finer
    wave_synth, flux_synth = _read_smooth('smooth.out')
    wavelength_obs, flux_obs = read_observations(wavelength, flux, start_synth, end_synth)
    f = interp1d(wave_synth, flux_synth, kind='cubic')
    flux_inter_synth = f(wavelength_obs)
    return wavelength_obs, flux_obs, flux_inter_synth


def plot_synthetic():
    """Function to plot synthetic spectrum
    """
    import seaborn
    import matplotlib.pyplot as plt
    x, y = _read_smooth(fname='smooth.out')
    plt.plot(x, y)
    plt.show()
    plt.close()
