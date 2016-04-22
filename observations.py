# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import pylab as pl
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from glob import glob
import os

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
    """Read the file that contains the line list then read the lines"""

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


def _read_raw_moog(fname='summary.out'):
    """Read the summary.out and return them

    :fname: From the summary_out
    :returns: flux
    """
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
    """Read the synthetic spectrum from the summary.out and return them
    :fname: From the summary_out
    :returns: wavelength and flux
    """
    with open(fname, 'r') as f:
         lines = f.readlines()[2:]

    wavelength = []
    flux = []
    for line in lines:
        line = filter(None, line.split())
        wavelength.append(float(line[0]))
        flux.append(float(line[1]))
    return wavelength, flux

def read_synth_intervals(fname):
    """Read the synthetic from all the intervals
    :fname: Line list that includes the line list files
    :returns: wavelength and flux
    """
    n_intervals, ranges, fout = read_linelist(fname)
    spec = []
    for i in range(n_intervals):
        x, y = _read_moog('results/%s.spec' % fout[i]) 
        spec.append(_read_moog('results/%s.spec' % fout[i]))

    wavelength = np.column_stack(spec)[0]
    flux = np.column_stack(spec)[1]
    return wavelength, flux

   
def read_observations(fname, start_synth, end_synth):
    """Read observed spectrum and return wavelength and flux"""
    import pyfits

    extension = ('.dat', '.txt', '.fits')
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

        wavelength_obs = wave[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]
        flux_obs = flux[np.where((wave >= float(start_synth)) & (wave <= float(end_synth)))]

    else: 
        print('Spectrum is not in acceptable format. Convert to ascii of fits.')
        wavelength_obs, flux_obs = (None, None)
    return wavelength_obs, flux_obs     

def observations(fname):
    """Read observed spectrum and return wavelength, flux and header"""
    import pyfits

    extension = ('.dat', '.txt', '.fits')
    if fname.endswith(extension):
        if fname[-4:]=='.dat' or fname[-4:]=='.txt':
            with open(fname, 'r') as f:
                lines = (line for line in f if not line[0].isalpha()) #skip header
                wave, flux = np.loadtxt(lines, unpack=True, usecols=(0, 1))
                header = None
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

    else: 
        print('Spectrum is not in acceptable format. Convert to ascii of fits.')
        wave, flux, header = (None, None, None)
    return wave, flux, header     


def interpol_synthetic(wavelength, flux, start_synth, end_synth):
    """Interpolation of the synthetic flux to the observed wavelength"""
    from scipy.interpolate import interp1d
    # The synthetic spectrum should be always finer
    wave_synth, flux_synth = _read_smooth('smooth.out')
    wavelength_obs, flux_obs = read_observations(wavelength, flux, start_synth, end_synth)
    f = interp1d(wave_synth, flux_synth, kind='cubic')
    flux_inter_synth = f(wavelength_obs)
    return wavelength_obs, flux_obs, flux_inter_synth


def plot_synth(fname): 
    """Function to plot synthetic spectrum
    """
    import matplotlib.pyplot as plt
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
    import matplotlib.pyplot as plt
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

    import matplotlib.pyplot as plt
    pl.plot(x, y)
    pl.plot(x_obs, y_obs)
    pl.show()
    pl.close()
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

def chi2(wave_obs, flux_obs, fname):
    """Some x2 statistics"""
    wave_synth, flux_synth = _read_moog(fname)
    #Interpolation of the synthetic flux to the observed wavelength
    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    f = interp1d(wave_synth, flux_synth, kind='cubic')
    int_flux = f(wave_obs)
    error = 1./200 
    chi = ((flux_obs - int_flux)/error)**2
    chi2 = np.sum(((flux_obs - int_flux)/error)**2)
    #which or the 2 is the correct format?
    chi = (flux_obs - int_flux)**2
    chi2 = np.sum((flux_obs - int_flux)**2)
    print('This is your chi2 value: '), chi2
    return

def instrumental(x,y,resolution):
    """Convolve synthetic spectrum with gaussian profile"""
    z = pyasl.instrBroadGaussFast(x, y, resolution, edgeHandling="firstlast")
    return x, z

def vmac():
    return

def vmic():
    return
