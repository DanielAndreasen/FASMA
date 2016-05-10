# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from astropy.io import fits

def save_synth_spec(x, y, fname='initial.spec'):
    '''Save synthetic spectrum of all intervals

    Input
    ----
    x : ndarray
      Wavelength
    y : ndarray
      Flux
    fname : str
      Filename of fits file

    Output
    -----
    fname fits file
    '''

    tbhdu = fits.BinTableHDU.from_columns([fits.Column(name='wavelength', format='D', array=x),
                                           fits.Column(name='flux', format='D', array=y)])
    tbhdu.writeto('results/%s' % fname, clobber=True)
    print('Synthetic spectrum saved: results/%s' % fname)


def broadening(x, y, vsini, vmac, resolution=None, epsilon=0.60):
    '''This function broadens the given data using velocity kernels,
    e.g. instrumental profile, vsini and vmac.
    Based on http://www.hs.uni-hamburg.de/DE/Ins/Per/Czesla/PyA/PyA/pyaslDoc/aslDoc/broadening.html
    Input
    ----
    x : ndarray
      wavelength
    y : ndarray
      flux
    resolution : float
      Instrumental resolution (lambda /delta lambda)
    vsini : float
      vsini in km/s
    vmac : float
      vmac in km/s

    Output
    -----
    y_broad : ndarray
      Broadened flux
    x : ndarray
      Same wavelength
    '''

    from PyAstronomy import pyasl
    from scipy.signal import fftconvolve
    from scipy.integrate import quad


    def instrumental_profile(x, y, resolution):
        '''
        Inputs
        -----
        x, y : The abscissa and ordinate of the data.
        sigma : The width (i.e., standard deviation) of the Gaussian profile
        used in the convolution.
        edgeHandling : None, "firstlast". Determines the way edges will be
        handled. If None, nothing will be done about it. If set to "firstlast",
        the spectrum will be extended by using the first and last value at the
        start or end. Note that this is not necessarily appropriate.
        The default is None.
        maxsig : The extent of the broadening kernel in terms of standrad
        deviations. By default, the Gaussian broadening kernel will be extended
        over the entire given spectrum, which can cause slow evaluation in the
        case of large spectra. A reasonable choice could, e.g., be five.

        Output
        -----
        y_inst : convolved flux
        '''

        #Deal with zero or None values seperately
        if (resolution is None) or (resolution == 0):
            y_inst = y
        else:
            y_inst = pyasl.instrBroadGaussFast(x, y, resolution, edgeHandling="firstlast", fullout=False, maxsig=None)
        return y_inst

    def vsini_broadening(x, y_inst, epsilon, vsini):
        '''
        Apply rotational broadening to a spectrum assuming a linear limb darkening
        law. The adopted limb darkening law is the linear one, parameterize by the
        linear limb darkening parameter: epsilon = 0.6.
        The effect of rotational broadening on the spectrum is
        wavelength dependent, because the Doppler shift depends
        on wavelength. This function neglects this dependence, which
        is weak if the wavelength range is not too large.
        Code from: http://www.phoebe-project.org/2.0/
        .. note:: numpy.convolve is used to carry out the convolution
              and "mode = same" is used. Therefore, the output
              will be of the same size as the input, but it
              will show edge effects.
        Input
        -----
        wvl : The wavelength
        flux : The flux
        epsilon : Linear limb-darkening coefficient (0-1).
        vsini : Projected rotational velocity in km/s.
        effWvl : The wavelength at which the broadening kernel is evaluated.
        If not specified, the mean wavelength of the input will be used.

        Output
        ------
        y_rot : convolved flux
        '''
        if vsini == 0:
            y_rot = y
        else:
            y_rot = pyasl.fastRotBroad(x, y, epsilon, vsini, effWvl=None)
        return y_rot

    def _vmacro_kernel(dlam, Ar, At, Zr, Zt):
        '''
        Macroturbulent velocity kernel.
        '''
        dlam[dlam == 0] = 1e-8
        if Zr != Zt:
            return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) * quad(lambda u: np.exp(-1/u**2),0,Zr/idlam)[0] + \
                          2*At*idlam/(np.sqrt(np.pi)*Zt**2) * quad(lambda u: np.exp(-1/u**2),0,Zt/idlam)[0])
                             for idlam in dlam])
        else:
            return np.array([(2*Ar*idlam/(np.sqrt(np.pi)*Zr**2) + 2*At*idlam/(np.sqrt(np.pi)*Zt**2))\
                           * quad(lambda u: np.exp(-1/u**2),0,Zr/idlam)[0]\
                             for idlam in dlam])

    def _broadening_macroturbulent(wave, flux, vmacro_rad, vmacro_tan=None,
                                   return_kernel=False):
        '''
        Apply macroturbulent broadening.
        The macroturbulent kernel is defined as in [Gray2005]:

        .. math::
            K_\mathrm{macro}(\Delta\lambda) = \frac{2A_R\Delta\lambda}{\sqrt{\pi}\zeta_R^2}\int_0^{\zeta_R/\Delta\lambda}e^{-1/u^2}du

             & + \frac{2A_T\Delta\lambda}{\sqrt{\pi}\zeta_T^2}\int_0^{\zeta_T/\Delta\lambda}e^{-1/u^2}du

        If :envvar:`vmacro_tan` is :envvar:`None`, then the value will be put equal
        to the radial component :envvar:`vmacro_rad`.

        Input
        -----
        :parameter wave: Wavelength of the spectrum
        :type wave: array
        :parameter flux: Flux of the spectrum
        :type flux: array
        :parameter vmacro_rad: macroturbulent broadening, radial component
        :type vmacro_rad: float
        :parameter vmacro_tan: macroturbulent broadening, tangential component
        :type vmacro_tan: float
        :parameter return_kernel: return kernel
        :type return_kernel: bool

        Output
        ------
        y_mac : broadened flux [, (wavelength, kernel)]
        '''

        if vmacro_tan is None:
            vmacro_tan = vmacro_rad

        if vmacro_rad == vmacro_tan == 0:
            return flux

        # Define central wavelength
        lambda0 = (wave[0] + wave[-1]) / 2.0

        vmac_rad = vmacro_rad/(299792458.*1e-3)*lambda0
        vmac_tan = vmacro_tan/(299792458.*1e-3)*lambda0

        # Make sure the wavelength range is equidistant before applying the
        # convolution
        delta_wave = np.diff(wave).min()
        range_wave = wave.ptp()
        n_wave = int(range_wave/delta_wave)+1
        wave_ = np.linspace(wave[0], wave[-1], n_wave)
        flux_ = np.interp(wave_, wave, flux)
        dwave = wave_[1]-wave_[0]
        n_kernel = int(5*max(vmac_rad, vmac_tan)/dwave)
        if n_kernel % 2 == 0:
            n_kernel += 1

        # The kernel might be of too low resolution, or the the wavelength range
        # might be too narrow. In both cases, raise an appropriate error
        if n_kernel == 0:
            raise ValueError(("Spectrum resolution too low for "
                          "macroturbulent broadening"))
        elif n_kernel > n_wave:
            raise ValueError(("Spectrum range too narrow for "
                          "macroturbulent broadening"))

        # Construct the broadening kernel
        wave_k = np.arange(n_kernel)*dwave
        wave_k -= wave_k[-1]/2.
        kernel = _vmacro_kernel(wave_k, 1.0, 1.0, vmac_rad, vmac_tan)
        kernel /= sum(kernel)

        flux_conv = fftconvolve(1-flux_, kernel, mode='same')
        # And interpolate the results back on to the original wavelength array,
        # taking care of even vs. odd-length kernels
        if n_kernel % 2 == 1:
            offset = 0.0
        else:
            offset = dwave / 2.0
        flux = np.interp(wave+offset, wave_, 1-flux_conv)

        # Return the results.
        if return_kernel:
            return flux, (wave_k, kernel)
        else:
            return flux

    # Instrumental broadening
    y_inst = instrumental_profile(x, y, resolution)
    # vsini broadening
    y_rot = vsini_broadening(x, y_inst, epsilon, vsini)
    # vmac broadening
    y_broad = _broadening_macroturbulent(x, y_rot, vmacro_rad=vmac,
                                         vmacro_tan=None, return_kernel=False)
    return x, y_broad


def _read_raw_moog(fname='summary.out'):
    '''Read the summary.out and return them

    :fname: From the summary_out
    :returns: flux
    '''
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
    '''Read the output of moog - synthetic spectra.
    Input
    -----
    fname: smooth.out

    Output
    ------
    wavelength
    flux
    '''

    wavelength, flux = np.loadtxt(fname, skiprows=2, usecols=(0, 1), unpack=True)
    return wavelength, flux

    return wavelength, flux


def read_wave(linelist):
    '''Read the wavelenth intervals of the line list'''

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
    '''Read the file that contains the line list then read the lines
    Input
    -----
    fname : file that contains the filenames of the linelist

    Output
    ------
    ranges : wavelength ranges of the linelist
    flines : filename of the linelist
    '''

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
    return ranges, flines


def read_moog_intervals(fname, options):
    '''Read all the synthetic spectra from all the intervals
    and save it in a spec (fits) file.
    Input:
    -----
    fname: Line list that includes the line list files

    Output:
    ------
    wavelength
    flux
    '''

    ranges, fout = read_linelist(fname)
    n_intervals = len(ranges)
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
    '''Interpolation of the synthetic flux to the observed wavelength.
    Input
    -----
    wave_obs : observed wavelength
    wave_synth : synthetic wavelength
    flux_synth : synthetic flux

    Output
    ------
    int_flux : Interpolated synthetic flux
    '''

    from scipy.interpolate import InterpolatedUnivariateSpline

    sl = InterpolatedUnivariateSpline(wave_synth, flux_synth, k=1)
    int_flux = sl(wave_obs)
    return int_flux
