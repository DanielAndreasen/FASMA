# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from astropy.io import fits as pyfits

def save_synth_spec(x, y, fname='initial.spec'):
    """Save synthetic spectrum of all intervals
    Input
    ----
    x : wavelength
    y : flux
    fname : filename of fits file

    Output
    -----
    fname fits file
    """

    tbhdu = pyfits.BinTableHDU.from_columns([pyfits.Column(name='wavelength', format='D', array=x), pyfits.Column(name='flux', format='D', array=y)])
    tbhdu.writeto('results/%s' % fname, clobber=True)
    print('Synthetic spectrum saved: results/%s' % fname)
    return


def broadening(x, y, resolution=None, vsini=0.0, epsilon=0.60, vmac=0.0):
    """This function broadens the given data using velocity kernels,
    e.g. instrumental profile, vsini and vmac.
    Based on http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/broadening.html
    Input
    ----
    x : wavelength
    y : flux
    resolution : instrumental resolution (lambda /delta lambda)
    vsini : in km/s
    vmac : in km/s

    Output
    -----
    y_broad : broadened flux
    x : same wavelength
    """

    from PyAstronomy import pyasl

    #instrumental broadening
    if (resolution is None) or (resolution == 0):
        y_inst = y
    else:
        y_inst = pyasl.instrBroadGaussFast(x, y, int(resolution), edgeHandling="firstlast", fullout=False, maxsig=None)

    #vsini broadening
    #Apply rotational broadening to a spectrum assuming a linear limb darkening
    #law. The adopted limb darkening law is the linear one, parameterize by the
    #linear limb darkening parameter: epsilon = 0.6.
    if vsini == 0:
        y_rot = y_inst
    else:
        y_rot = pyasl.fastRotBroad(x, y_inst, epsilon, vsini, effWvl=None)

    #vmac broadening
    if vmac == 0:
        y_broad = y_rot
    else:
       #I donow
       y_broad = y_rot
    return x, y_broad


def _read_raw_moog(fname='summary.out'):
    """Read the summary.out and return them

    :fname: From the summary_out
    :returns: flux
    """
    import itertools

    with open('summary.out', 'r') as f:
        lines = f.readlines()[:]

    start_wave, end_wave, step, flux_step = lines[2].split()
    flux_raw = []
    for line in lines[3:]:
       line = line.split()
       flux_raw.append(line)

    flux = np.array(list(itertools.chain.from_iterable(flux_raw)))
    flux = flux.astype(np.float)
    flux = 1.0-flux
    w0, dw, n = float(start_wave), float(step), len(flux)
    w = w0 + dw * len(flux)
    wavelength = np.linspace(w0, w, n, endpoint=False)
    return wavelength, flux


def _read_moog(fname='smooth.out'):
    """Read the output of moog - synthetic spectra.
    Input
    -----
    fname: smooth.out

    Output
    ------
    wavelength
    flux
    """

    wavelength, flux = np.loadtxt(fname, skiprows=2, usecols=(0, 1), unpack=True)
    return wavelength, flux

    return wavelength, flux


def read_wave(linelist):
    """Read the wavelenth intervals of the line list"""

    with open(linelist, 'r') as f:

        lines = f.readlines()
    first_line = lines[0].split()

    if len(first_line) == 1:
        start_wave = first_line[0].split('-')[0]
        end_wave = first_line[0].split('-')[1]
    else:
        start_wave = first_line[0]
        end_wave = lines[-1].split()[0]
    return start_wave, end_wave


def read_linelist(fname):
    """Read the file that contains the line list then read the lines
    Input
    -----
    fname : file that contains the filenames of the linelist

    Output
    ------
    n_intervals : Number of intervals in the linelist
    ranges : wavelength ranges of the linelist
    flines : filename of the linelist
    """

    with open('linelist/%s' % fname, 'r') as f:

        lines = f.readlines()

    n_intervals = len(lines)
    ranges = []
    flines = []
    for line in lines:
        line = line.split()
        # Check if linelist files are inside the directory, if not break
        if not os.path.isfile('linelist/%s' % line[0]):
            raise IOError('The linelist is not in the linelist directory!')
        flines.append(line[0])

        with open('linelist/%s' % line[0], 'r') as f:

            lines = f.readlines()
        first_line = lines[0].split()

        if len(first_line) == 1:
            start_wave = first_line[0].split('-')[0]
            end_wave = first_line[0].split('-')[1]
            r = (float(start_wave), float(end_wave))
            ranges.append(r)
        else:
            start_wave = first_line[0]
            end_wave = lines[-1].split()[0]
            r = (float(start_wave), float(end_wave))
            ranges.append(r)
    return n_intervals, ranges, flines


def read_moog_intervals(fname, options):
    """Read all the synthetic spectra from all the intervals
    and save it in a spec (fits) file.
    Input:
    -----
    fname: Line list that includes the line list files

    Output:
    ------
    wavelength
    flux
    """

    n_intervals, ranges, fout = read_linelist(fname)
    spec = []
    for i in range(n_intervals):
        x, y = _read_moog('results/%s.spec' % fout[i])
        spec.append(_read_moog('results/%s.spec' % fout[i]))

    w = np.column_stack(spec)[0]
    f = np.column_stack(spec)[1]

    #Add broadening
    wavelength, flux = broadening(w, f, resolution=options['resolution'],
    vsini=options['vsini'], epsilon=options['limb'], vmac=options['vmac'])
    return wavelength, flux


def interpol_synthetic(wave_obs, wave_synth, flux_synth):
    """Interpolation of the synthetic flux to the observed wavelength.
    Input
    -----
    wave_obs : observed wavelength
    wave_synth : synthetic wavelength
    flux_synth : synthetic flux

    Output
    ------
    int_flux : Interpolated synthetic flux
    """

    from scipy.interpolate import InterpolatedUnivariateSpline

    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux
