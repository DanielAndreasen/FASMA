# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from glob import glob
import os
from astropy.io import fits as pyfits


def read_observations(fname, start_synth, end_synth):
    """Read observed spectrum and return wavelength and flux"""

    extension = ('.dat', '.txt', '.spec', '.fits')
    if fname.endswith(extension):
        if fname[-4:]=='.dat' or fname[-4:]=='.txt':
            with open(fname, 'r') as f:
                lines = (line for line in f if not line[0].isalpha()) #skip header
                wave, flux = np.loadtxt(lines, unpack=True, usecols=(0, 1))

        elif fname[-5:]=='.fits':
            hdulist = pyfits.open(fname)
            header = hdulist[0].header
            flux = hdulist[0].data #flux data in the primary
            flux = np.array(flux, dtype=np.float64)
            start_wave = header['CRVAL1'] #initial wavelenght
            #step = header['CD1_1'] #step in wavelenght
            step = header['CDELT1'] #increment per pixel
            w0, dw, n = start_wave, step, len(flux)
            w = start_wave + step * n
            wave = np.linspace(w0, w, n, endpoint=False)

        elif fname[-5:]=='.spec':
            hdulist = pyfits.open(fname)
            x = hdulist[1].data
            flux = x['flux'] #flux data in the primary
            wave = x['wavelength']

        wavelength_obs = wave[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
        flux_obs = flux[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]

    else:
        print('Spectrum is not in acceptable format. Convert to ascii of fits.')
        wavelength_obs, flux_obs = (None, None)
    return wavelength_obs, flux_obs


def read_obs_intervals(obs_fname, N, r):
    """Read only the spectral chunks from the observed spectrum"""
    spec = []
    for i in range(N):
        spec.append(read_observations(obs_fname, start_synth=r[i][0], end_synth=r[i][1]))

    x_obs = np.column_stack(spec)[0]
    y_obs = np.column_stack(spec)[1]
    return x_obs, y_obs


def plot_synth(fname):
    """Function to plot synthetic spectrum
    """

    if fname=='smooth.out':
        x, y = _read_moog(fname='smooth.out')
        pl.plot(x, y)
        pl.show()
        pl.close()
    elif fname=='summary.out':
        x, y = _read_raw_moog('summary.out')
        #z = pyasl.instrBroadGaussFast(x, y, 50000, edgeHandling="firstlast")
        pl.plot(x, y)
        #pl.plot(x, z)
        pl.show()
        pl.close()
    else:
        if os.path.isfile(fname):
            x, y = _read_moog(fname)
            pl.plot(x, y)
            pl.show()
            pl.close()
        else:
            print('Synthetic spectrum does not exist.')
    return


def plot_synth_obs(x_obs, y_obs, fname):
    """Function to plot synthetic spectrum
    """

    if fname=='smooth.out':
        x, y = _read_moog(fname='smooth.out')
        pl.plot(x, y)
        pl.plot(x_obs, y_obs)
        pl.show()
        pl.close()
    elif fname=='summary.out':
        x, y = _read_raw_moog('summary.out')
        #z = pyasl.instrBroadGaussFast(x, y, 50000, edgeHandling="firstlast")
        pl.plot(x, y)
        pl.plot(x_obs, y_obs)
        #pl.plot(x, z)
        pl.show()
        pl.close()
    else:
        if os.path.isfile(fname):
            x, y = _read_moog(fname)
            pl.plot(x, y, '-')
            pl.plot(x_obs, y_obs, '--')
            pl.show()
            pl.close()
        else:
            print('Synthetic spectrum does not exist.')
    return


def plot(x_obs, y_obs, x, y):
    """Function to plot synthetic spectrum
    """

    if (x_obs is None) and (x is None):
        pass
    if x_obs is None:
        pl.plot(x, y)
        pl.show()
    else:
        pl.plot(x, y)
        pl.plot(x_obs, y_obs)
        pl.show()
    return


def interpol_synthetic(wave_obs, wave_synth, flux_synth):
    """Interpolation of the synthetic flux to the observed wavelength"""
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux


def normalize(x, y, num=20, k=1):
    index_max = np.sort(np.argsort(y)[-num:])
    f_obs_max = y[index_max]
    w_obs_max = x[index_max]
    sl = InterpolatedUnivariateSpline(w_obs_max, f_obs_max, k=1)
    continuum = sl(x)
    norm = y/continuum
    return x, norm


def chi2(wave_obs, flux_obs, wave_synth, flux_synth):
    """Some x2 statistics"""

    #Interpolation of the synthetic flux to the observed wavelength
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    #f = interp1d(wave_synth, flux_synth, kind='cubic')
    #error = 1./200
    #chi2 = np.sum(((flux_obs - int_flux)/error)**2)
    #which or the 2 is the correct format?
    chi = ((flux_obs - int_flux)**2)
    chi2 = np.sum(chi)
    print('This is your chi2 value: '), chi2
    return
