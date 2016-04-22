#!/usr/bin/env python
# -*- coding: utf8 -*-
import numpy as np
from astropy.io import fits as pyfits
import argparse
import pylab as pl
from scipy import stats
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from observations import observations
import logging
import os
from matplotlib import pyplot as pl

def _options(options=None):
    """Reads the options inside the config file"""
    defaults = {'output': False,
                'chunks': False,
                'iterations':1,
                'cosmic': False, 
                'clip':4.0,
                'plot': False
                }
    if not options:
        return defaults
    else:
        for option in options.split(','):
            if ':' in option:
                option = option.split(':')
                defaults[option[0]] = option[1]
            else:
                if option in ['chunks', 'plot', 'cosmic']:
                    defaults[option] = False if defaults[option] else True
        defaults['clip'] = float(defaults['clip'])
        defaults['iterations'] = int(defaults['iterations'])
        defaults['chunks'] = int(defaults['chunks'])
        return defaults


def remove_cosmic(filename, clip=4.0, num_chunks=50, plot=False):
    """Remove cosmic rays with sigma clipping."""

    def mad(data, axis=None):
        return np.median(np.absolute(data - np.median(data, axis)), axis)

    def clean(chunks, clip):
        """Divide into chunks and remove outliers mad*clip from the median flux""" 
        wave = chunks[0]
        flux = chunks[1]
        med = np.median(flux)
        sigma = mad(flux)
        n = len(flux)
        bad_points = flux[np.where(flux> (med + (sigma*clip)))]
        bad_wave = wave[np.where(flux> (med + (sigma*clip)))]
        fluxout = np.zeros(n)   
        for i in range(n):
            if flux[i] < (med + (sigma*clip)):
	        fluxout[i] = flux[i]   
            if flux[i] > (med + (sigma*clip)):
                fluxout[i-2:i+2] = med #this should be done better or iteratively

        return wave, fluxout, bad_wave, bad_points

    wave_obs, flux_obs, z = observations(filename)
    #Ignore zero or negative flux points
    w0 = wave_obs[np.where(flux_obs <= 0)]
    f0 = flux_obs[np.where(flux_obs <= 0)]
    index = np.where(flux_obs <= 0)[0]

    wave = wave_obs[np.where(flux_obs > 0)]
    flux = flux_obs[np.where(flux_obs > 0)]

    #Divide spectrum in chunks
    spectrum = np.array([wave, flux])
    chunks = np.array_split(spectrum, num_chunks, axis=1)
    clean_flux = np.zeros(len(chunks))
    new_flux = []
    new_wave = []
    new_xbad = []
    new_ybad = []
    for chunk in chunks:
        wave, clean_flux, x_bad, y_bad = clean(chunk, clip)
        for f in clean_flux:
            new_flux.append(float(f))
        for w in wave:
            new_wave.append(float(w))
        for xb in x_bad:
            new_xbad.append(float(xb))
        for yb in y_bad:
            new_ybad.append(float(yb))

    for x in reversed(w0): 
        new_wave.insert(index[0],x)
    for x in reversed(f0): 
        new_flux.insert(index[0],x)

    flux  = np.array(new_flux, dtype=np.float64)
    wave  = np.array(new_wave, dtype=np.float64)
    x_bad = np.array(new_xbad, dtype=np.float64)
    y_bad = np.array(new_ybad, dtype=np.float64)

    if plot: 
        pl.plot(wave, flux)
        pl.plot(x_bad,y_bad, 'o')
        pl.show()

    return wave, flux

def normalize(wave_obs, flux_obs, num_chunks, plot=False):
    """Normalise function. Make a polynomial fit from the maximum points of each segment"""

    #Ignore zero or negative flux points
    w0 = wave_obs[np.where(flux_obs <= 0)]
    f0 = flux_obs[np.where(flux_obs <= 0)]
    index = np.where(flux_obs <= 0)[0]

    wave = wave_obs[np.where(flux_obs > 0)]
    flux = flux_obs[np.where(flux_obs > 0)]

    spectrum = np.array([wave, flux])
    chunks = np.array_split(spectrum, num_chunks, axis=1)

    clean_flux = np.zeros(len(chunks))
    max_flux = []
    max_wave = []
    for chunk in chunks:
        wave_c = chunk[0]
        flux_c = chunk[1]
        index_max = np.sort(np.argsort(flux_c)[-8:]) # this can be done better
        f_max = flux_c[index_max]
        w_max = wave_c[index_max]
        for x in f_max:
            max_flux.append(float(x))
        for y in w_max:
             max_wave.append(float(y))
        
    f = np.array(max_flux, dtype=np.float64)
    w = np.array(max_wave, dtype=np.float64)

    z = np.polyfit(max_wave, max_flux, 5)
    p = np.poly1d(z)
    #Taking care of the limits
    f[0] = p(wave)[0]
    w[0] = wave[0]
    f[-1] = p(wave)[-1]
    w[-1] = wave[-1]

    sl = InterpolatedUnivariateSpline(w, f, k=1)
    new_flux = flux/sl(wave)

    #Add here the negative values
    wave = np.insert(wave, index[0], w0)
    new_flux = np.insert(new_flux, index[0], f0)

    if plot:
        pl.plot(wave_obs, flux_obs)
        pl.plot(wave_obs, sl(wave_obs))
        pl.plot(max_wave, max_flux, 'o')
        pl.show()
    
        pl.plot(wave_obs, flux_obs/sl(wave_obs), label='spline')
        pl.legend()
        pl.show()

    return wave, new_flux


def normdriver(starLines='StarMe_norm.cfg', overwrite=None):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe.cfg)
    overwrite   -   Overwrite the results.csv file (default: False)

    """
    try:  # Cleaning from previous runs
        os.remove('captain.log')
    except OSError:
        pass
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler('captain.log')
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalnum():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Spectrum: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            if len(line) not in [1, 2]:
                logger.error('Could not process information for this spectrum: %s' % line)
                continue

            # Check if the spectrum exists, if not log it and pass to next
            if not os.path.isfile('spectra/%s' % line[0]):
                logger.error('Error: %s not found.' % line[0])
                continue

            if len(line) == 1:
                if os.path.isfile('spectra/%s' % line[0]): 
                    spectrum = 'spectra/%s' % line[0]
                elif os.path.isfile(line[0]):
                    spectrum = line[0]
                else: 
                    print('Error: %s not found.' % line[0])

                options = _options()
                #no need to clean for cosmics
                x_obs, y_obs, header = observations(spectrum)
                n = (x_obs[-1]-x_obs[0])/20. #divided in chuncks of 20A
                n = int(n)
                x_n, y_n = normalize(x_obs, y_obs, n)

            elif len(line) == 2:
                if os.path.isfile('spectra/%s' % line[0]): 
                    spectrum = 'spectra/%s' % line[0]
                elif os.path.isfile(line[0]):
                    spectrum = line[0]
                else: 
                    print('Error: %s not found.' % line[0])

                options = _options(line[-1])

                if options['cosmic']: #clean for cosmic rays
                    if options['plot']:
                        x, y = remove_cosmic(spectrum, clip=options['clip'], num_chunks=50, plot=True)
                    else: 
                        x, y = remove_cosmic(spectrum, clip=options['clip'], num_chunks=50, plot=False)
            
                    if options['chunks']:
                        x_n, y_n = normalize(x, y, options['chunks'], plot=True)
                    else: 
                        x_obs, y_obs, header = observations(spectrum)
                        n = (x_obs[-1]-x_obs[0])/20. #divided in chuncks of 20A
                        n = int(n)
                        x_n, y_n = normalize(x, y, n, plot=True)

                else: #no need to clean for cosmics
                    x_obs, y_obs, header = observations(spectrum)
                    if options['plot']:
                        if options['chunks']:
                            print(options)
                            x_n, y_n = normalize(x_obs, y_obs, options['chunks'], plot=True)
                        else: 
                            n = (x[-1]-x[0])/20. #divided in chuncks of 20A
                            n = int(n)
                            x_n, y_n = normalize(x_obs, y_obs, n, plot=True)
                    else:
                        if options['chunks']:
                            x_n, y_n = normalize(x_obs, y_obs, options['chunks'], plot=False)
                        else: 
                            n = (x[-1]-x[0])/20. #divided in chuncks of 20A
                            n = int(n)
                            x_n, y_n = normalize(x_obs, y_obs, n, plot=False)
                prihdr = pyfits.Header()
                prihdr["NAXIS1"] = len(x_n)
                prihdr["CDELT1"] = x_n[1]-x_n[0]
                prihdr["CRVAL1"] = x_n[0]
                pyfits.writeto(spectrum+'_norm.fits', y_n, prihdr, clobber=True)
                #np.savetxt(spectrum+'_norm.txt', zip(x_n, y_n), delimiter='\t')
                print('The normalized spectrum is here: %s_norm.fits' % spectrum)
    return x_n, y_n

if __name__ == '__main__':
    x_n, y_n = normdriver()

