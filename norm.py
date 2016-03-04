#!/usr/bin/python
import numpy as np
from astropy.io import fits as pyfits
import argparse
import pylab as pl
import scipy
from scipy import stats
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

parser = argparse.ArgumentParser(description='Normalise')
parser.add_argument('input', help='Input file')
args = parser.parse_args()

def read_fits(filename):
    hdulist = pyfits.open(filename)
    header = hdulist[0].header #read header
    flux = hdulist[0].data #flux data in the primary
    flux = np.array(flux, dtype=np.float64)
    start_wave = header['CRVAL1'] #initial wavelenght
    step = header['CD1_1'] #step in wavelenght
    w0, dw, n = start_wave, step, len(flux)
    w = start_wave + step * n
    wave = np.linspace(w0, w, n, endpoint=False)
    return wave, flux, header

def rv_correction(filename, velocity):
    c = 299792.458 #km s-1 
    hdulist = pyfits.open(filename)
    header = hdulist[0].header #read header
    flux = hdulist[0].data #flux data in the primary
    flux = np.array(flux, dtype=np.float64)
    start_wave = header['CRVAL1'] #initial wavelenght
    step = header['CD1_1'] #step in wavelenght
    # New starting wavelenght and step
    relat_wave = start_wave * np.sqrt((1. - velocity/c)/(1. + velocity/c)) 
    relat_step = step * np.sqrt((1. - velocity/c)/(1. + velocity/c)) 
    #Change the header with the new values 
    header['CRVAL1'] = relat_wave
    header['CDELT1'] = relat_step
    header['CD1_1'] = relat_step
    #Save results
    pyfits.writeto(filename[:-4]+"_rv.fits", flux, header)
    hdulist.close()
    return

def plot(filename):
    wave, flux, z = read_fits(filename)
    pl.plot(wave, flux)
    pl.show()
    return

def mad(data, axis=None):
    return np.median(np.absolute(data - np.median(data, axis)), axis)

def remove_cosmics(chunks, clip):
    wave = chunks[0]
    flux = chunks[1]
    med = np.median(flux)
    sigma = mad(flux)
    n = len(flux)
    bad_points = flux[np.where(flux> (med + (sigma*clip)))]
    x = wave[np.where(flux> (med + (sigma*clip)))]
    fluxout = np.zeros(n)   
    for i in range(n):
        if flux[i] < (med + (sigma*clip)):
	   fluxout[i] = flux[i]   
        if flux[i] > (med + (sigma*clip)):
           fluxout[i-2:i+6] = med #this should be done better
    pl.plot(x, bad_points, 'o')
    return wave, fluxout

def fits_chunk_division(filename, num_chunks):
    """Divided the fits file into a defined amound of chunks."""
    wave, flux, z= read_fits(filename)
    n_flux = len(flux)
    # remove negative values for flux
    flux = scipy.stats.threshold(flux, threshmin=0.0, threshmax='None', newval=0.01)
    #Divide spectrum in chunks
    spectrum = np.array([wave, flux])
    chunks = np.array_split(spectrum, num_chunks, axis=1)
    clean_flux = np.zeros(len(chunks))
    new_flux = []
    new_wave = []
    for chunk in chunks:
        wave, clean_flux = remove_cosmics(chunk, clip=5.0)
        for sub_flux in clean_flux:
            new_flux.append(float(sub_flux))
        for sub_flux in wave:
            new_wave.append(float(sub_flux))

    flux = np.array(new_flux, dtype=np.float64)
    wave = np.array(new_wave, dtype=np.float64)
    return wave, flux

def norm(wave, flux, num_chunks):
    """Normalise after cosmic removal"""
    #make a polynomial fit from the maximum points of each segment
    spectrum = np.array([wave, flux])
    chunks = np.array_split(spectrum, num_chunks, axis=1)
    clean_flux = np.zeros(len(chunks))
    max_flux = []
    max_wave = []
    for chunk in chunks:
        wave_c = chunk[0]
        flux_c = chunk[1]
        index_max = np.sort(np.argsort(flux_c)[-8:])
        f_max = flux_c[index_max]
        w_max = wave_c[index_max]
        pl.plot(w_max, f_max, 'o')
        pl.plot(wave_c, flux_c)
        for x in f_max:
            max_flux.append(float(x))
        for y in w_max:
             max_wave.append(float(y))
        
    f = np.array(max_flux, dtype=np.float64)
    w = np.array(max_wave, dtype=np.float64)
    f = np.insert(f, 0, np.max(flux[0:10]))
    w = np.insert(w, 0, wave[0])
    f = np.insert(f, -1, np.max(flux[-10:-1]))
    w = np.insert(w, -1, wave[-1])
    z = np.polyfit(max_wave, max_flux, 5)
    p = np.poly1d(z)
    pl.plot(wave, flux)
    pl.plot(wave, p(wave), 'b')
    sl = scipy.interpolate.InterpolatedUnivariateSpline(w, f, k=1)
    pl.plot(wave, sl(wave))
    pl.plot(max_wave, max_flux, 'o')
    pl.show()
    
    pl.plot(wave, flux/p(wave))
    pl.plot(wave, flux/sl(wave), label='spline')
    pl.legend()
    #plot('SUN/SUN_norm.fits')
    pl.show()
    new_flux = flux/sl(wave)
    return wave, new_flux

#Main
if __name__ == '__main__':

    x, y, h = read_fits(args.input)
    x_c, y_c = fits_chunk_division(args.input, num_chunks=5)
    pl.plot(x,y, x_c,y_c)
    pl.show()
    w, f, header = read_fits(args.input)
    x, y = fits_chunk_division(args.input, num_chunks=5)
    x_new, y_new = norm(x, y, num_chunks=40)
    pyfits.writeto(args.input[:-5]+"_norm.fits", y_new, header)

