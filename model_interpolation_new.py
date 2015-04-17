#!/usr/bin/python
from __future__ import division
import numpy as np
from scipy.integrate import simps
from scipy import integrate
import scipy as sp
import scipy
from scipy.interpolate import interp1d
import gzip
import matplotlib.pyplot as plt


"""
Following the concept of Sz. Mezeros and C. Allende Prieto: For each set
of parameters, we identified the 8 immediate neighbors with higher and
lower values for each parameter in the grid, calculated by numerical
integration the Rosseland optical depth for each, re-sampled all the
thermodynamical quantities in the atmosphere (temperature, gas pressure,
and electron density) on a common optical depth scale for all models
by linear interpolation, and then interpolated, linearly, all the
thermodynamical quantities to the parameters (Teff , log g, and [Fe/H])
of the target model. Other quantities included in the models (Rosseland
opacities, radiative pressure, etc.) were also interpolated in the same
way.
"""


def _unpack_model(fname):
    """Unpack the compressed model and store it in a temporary file

    :fname: File name of the compressed atmosphere model
    :returns: String of the uncompressed atmosphere model
    """
    f = gzip.open(fname)
    return f.readlines()


def _read_header(fname):
    """Read the header of the model atmosphere

    :fname: file name of the model atmosphere
    :returns: Teff and number of layers

    """
    # Get information from the header
    teff, num_layers = None, None
    # with open(fname) as lines:
    for line in fname:
        vline = line.split()
        param = vline[0]
        if param == 'TEFF':
            teff = float(vline[1])
            logg = float(vline[3])
        elif param == "READ":
            num_layers = int(vline[2])
        # Exit the loop when we have what we need
        if teff and num_layers:
            return teff, num_layers


def save_model(model, type='kurucz95', fout='out.atm', vt=1.2):
    """Save the model atmosphere in the right format

    :model: The interpolated model atmosphere
    :type: Type of model atmosphere (onyl kurucz95 at the moment)
    :fout: Which place to save to
    """
    if type == 'kurucz95':
        header = 'KURUCZ\n'\
                 'Teff=%i   log g=%.2f   [Fe/H]=%.2f    vt=%.3e\n'\
                 'NTAU        %i' % (5777, 4.44, -0.14, 2.4e5, 72)
    elif type.lower() == 'marcz':  # How to spell this?
        raise NotImplementedError('Patience is the key. Wait a bit more for %s\
                                   models to be implemented.' % type)
    else:
        raise NameError('Could not find %s models' % type)

    footer = '    %.3e\n'\
             'NATOMS     1  %.2f\n'\
             '      26.0   %.2f\n'\
             'NMOL      19\n'\
             '      606.0    106.0    607.0    608.0    107.0    108.0    112.0  707.0\n'\
             '       708.0    808.0     12.1  60808.0  10108.0    101.0     6.1    7.1\n'\
             '         8.1    822.0     22.1' % (vt*1e5, -0.2, 7.47-0.2)

    np.savetxt(fout, model.T, header=header, footer=footer, comments='',
               delimiter=' ',
               fmt=('%15.8E', '%8.1f', '%.3E', '%.3E', '%.3E', '%.3E', '%.3E'))


def tauross_scale(abross, rhox, num_layers):
    """Build the tau-ross scale

    :abross: absorption
    :rhox: density
    :num_layers: Number of layers in the model atmosphere
    :returns: the new tau-ross scale
    """

    tauross = sp.integrate.cumtrapz(rhox * abross, initial=rhox[0] * abross[0])

    # tauross = np.zeros(num_layers)
    # # This is supposed to be the first element
    # tauross[0] = abross[0] * rhox[0]
    # for i in range(2, num_layers+1):
    #     tauross[i-1] = sp.integrate.simps(rhox[0:i], abross[0:i],
    #                                         even='last')
    #     tauross[i-1] = np.trapz(rhox[0:i], abross[0:i])
    return tauross


def read_model(fname):
    """Read the model, return all the columns and tauross"""

    data = _unpack_model(fname)
    teff, num_layers = _read_header(data)

    # This are the thermodynamical quantities.
    model = np.genfromtxt(data, skiprows=23, skip_footer=2,
                          usecols=(0, 1, 2, 3, 4, 5, 6),
                          names=['RHOX', 'T', 'P', 'XNE', 'ABROSS',
                                 'ACCRAD', 'VTURB'])

    model_rhox = model['RHOX']
    model_t = model['T']
    model_p = model['P']
    model_xne = model['XNE']
    model_abross = model['ABROSS']
    model_accrad = model['ACCRAD']
    # TODO: We don't need this one. manual p. 17
    model_vturb = model['VTURB']

    tauross = tauross_scale(model_abross, model_rhox, num_layers)
    return (model_rhox, model_t, model_p, model_xne, model_abross,
            model_accrad, model_vturb, tauross)


def interp_model(tauross, model, tauross_new):
    """Interpolate a physical quantity from the model from the tauross scale to
    1 value of the new tauross scale

    :tauross: Old tauross scale
    :model: Column in the atmospheric model to be interpolated
    :tauross_new: New tauross scale
    :returns: The interpolated atmospheric model to the new scale
    """
    # Extra key-words speed up the function with a factor of 10!
    f = interp1d(tauross, model, assume_sorted=True, copy=False)
    # f = InterpolatedUnivariateSpline(tauross,model,k=3)
    return f(tauross_new)


def interpolator(mnames, teff, logg, feh, out='out.atm'):
    """The function to call from the main program (pymoog.py or main.py)

    :mnames: As generated from _get_models
    :out: The interpolated model saved in this file
    """

    teff, tefflow, teffhigh = teff[0], teff[1][0], teff[1][1]
    logg, logglow, logghigh = logg[0], logg[1][0], logg[1][1]
    feh, fehlow, fehhigh = feh[0], feh[1][0], feh[1][1]

    models = []
    opmin, opmax = [], []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm)
        opmin.append(np.amin(tatm[-1]))
        opmax.append(np.amax(tatm[-1]))

    coord = []
    coord.append([(teff-tefflow)/(teffhigh-tefflow)])
    coord.append([(logg-logglow)/(logghigh-logglow)])
    coord.append([(feh-fehlow)/(fehhigh-fehlow)])

    newdeck = np.zeros_like(models[0])
    print models[0][7]
    newdeck[:, 7] = models[0][:, 7]

    print newdeck
    print newdeck.shape









    return 0, 0, 0
