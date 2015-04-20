#!/usr/bin/python
from __future__ import division
import numpy as np
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from scipy import ndimage
import gzip


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


def _read_header(lines):
    """Read the header of the model atmosphere

    :lines: The lines from the atmosphere model (already a string at this
            point)
    :returns: Teff and number of layers

    """
    # Get information from the header
    teff, num_layers = None, None
    # with open(fname) as lines:
    for line in lines:
        vline = line.split()
        param = vline[0]
        if param == 'TEFF':
            teff = float(vline[1])
            # logg = float(vline[3])
        elif param == "READ":
            num_layers = int(vline[2])
        # Exit the loop when we have what we need
        if teff and num_layers:
            return teff, num_layers


def tauross_scale(abross, rhox, num_layers):
    """Build the tau-ross scale

    :abross: absorption
    :rhox: density
    :num_layers: Number of layers in the model atmosphere
    :returns: the new tau-ross scale
    """
    # x8 times faster than Jo Bovy's solution with int_newton_cotes
    tauross = cumtrapz(rhox * abross, initial=rhox[0] * abross[0])
    return tauross


def read_model(fname):
    """Read the model atmosphere

    :fname: The gz file of the model atmosphere
    :returns: The columns and tauross in a tuple

    """
    data = _unpack_model(fname)
    teff, num_layers = _read_header(data)
    # These are the thermodynamical quantities
    model = np.genfromtxt(data, skiprows=23, skip_footer=2,
                          usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8),
                          names=['RHOX', 'T', 'P', 'XNE', 'ABROSS',
                                 'ACCRAD', 'VTURB', 'tab1', 'tab2'])

    tauross = tauross_scale(model['ABROSS'], model['RHOX'], num_layers)
    return (model, tauross)


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
    return f(tauross_new)


def interpolator(mnames, teff, logg, feh, out='out.atm'):
    """The function to call from the main program (pymoog.py or main.py)

    :mnames: As generated from _get_models
    :teff: Requested Effective Temperature and the two closest models in
           the gridpoints
    :logg: Requested Surface gravity and the two closest models
    :feh: Requested metallicity and the two closest models
    :out: The interpolated model saved in this file
    """

    # TODO: This is super ugly way of doing stuff!
    teff, tefflow, teffhigh = teff[0], teff[1][0], teff[1][1]
    logg, logglow, logghigh = logg[0], logg[1][0], logg[1][1]
    feh, fehlow, fehhigh = feh[0], feh[1][0], feh[1][1]

    interGridShape = (2, 2, 2)

    # Reading the models and defining the opacity intervals
    models = []
    opmin, opmax = [], []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(np.array(tatm[0].tolist()))
        opmin.append(np.amin(tatm[1]))
        opmax.append(np.amax(tatm[1]))

    # Define the grid coordinates for the interpolation
    c1 = [(teff-tefflow)/(teffhigh-tefflow)]
    c2 = [(logg-logglow)/(logghigh-logglow)]
    c3 = [(feh-fehlow)/(fehhigh-fehlow)]
    coord = [c1, c2, c3]

    # Interpolate the models using the Force
    # Look at Jobovy code.
    newdeck = np.zeros_like(models[0])
    newdeck[:, 7:] = models[0][:, 7:]
    for layer in range(newdeck.shape[0]):
        for column in range(newdeck.shape[1]):
            tlayer = np.zeros(interGridShape)
            cntr = 0
            while cntr < np.prod(interGridShape):
                idx = np.unravel_index(cntr, interGridShape)
                tlayer[idx] = models[cntr][layer, column]
                cntr += 1
            newdeck[layer, column] =\
                ndimage.interpolation.map_coordinates(tlayer, coord, order=1)

    # ADD OUTPUT TO FILE
    return newdeck


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
