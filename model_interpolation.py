#!/usr/bin/python
from __future__ import division
import numpy as np
from operator import mul
from scipy.ndimage import _nd_image
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


def map_coordinates(input, coordinates, output=None, order=1,
                    mode='constant', cval=0.0):
    """
    This is a copy of the original code, however all checks are removed.
    That is okay, since we know what we are doing, and want a speedup.
    The code is 3x faster than standard scipy!

    See documentation for
        scipy.ndimage.map_coordinates
    """

    mode_dir = {'nearest': 0,
                'wrap': 1,
                'reflect': 2,
                'mirror': 3,
                'constant': 4}

    output_shape = coordinates.shape[1:]
    mode = mode_dir[mode]
    output = np.zeros(output_shape, dtype=input.dtype.name)
    return_value = output
    _nd_image.geometric_transform(input, None, coordinates, None, None,
                                  output, order, mode, cval, None, None)
    return return_value


def _unpack_model(fname):
    """Unpack the compressed model and store it in a temporary file

    :fname: File name of the compressed atmosphere model
    :returns: String of the uncompressed atmosphere model
    """
    f = gzip.open(fname, compresslevel=1)
    return f.readlines()


def _tupleset(t, i, value):
    """
    This is used in tauross_scale
    """
    l = list(t)
    l[i] = value
    return tuple(l)


def tauross_scale(abross, rhox, dx=1, axis=-1):
    """Build the tau-ross scale

    Note that this is the source of scipy.integrate.cumtrapz
    in a more raw format to speed things up.

    :abross: absorption
    :rhox: density
    :returns: the new tau-ross scale
    """

    y = rhox * abross
    d = dx
    initial = rhox[0] * abross[0]

    nd = len(y.shape)
    slice1 = _tupleset((slice(None),)*nd, axis, slice(1, None))
    slice2 = _tupleset((slice(None),)*nd, axis, slice(None, -1))
    res = np.add.accumulate(d * (y[slice1] + y[slice2]) / 2.0, axis)

    shape = list(res.shape)
    shape[axis] = 1
    res = np.concatenate([np.ones(shape, dtype=res.dtype) * initial, res],
                         axis=axis)

    # Original piece of code below. More elegant, but slower
    # tauross = cumtrapz(rhox * abross, initial=rhox[0] * abross[0])
    # return tauross
    return res


def _loadtxt(lines):
    """
    Super efficient 'loadtxt' that works exactly for our data.
    25x speed-up compared to fully np.loadtxt (which does a lot of checking)

    :return: numpy array of atmospheric quantities
    """
    row, col = len(lines), len(lines[0].split())
    data = np.empty((row, col))
    for rowi, line in enumerate(lines):
        data[rowi, :] = line.split()
    return data


def read_model(fname):
    """Read the model atmosphere

    :fname: The gz file of the model atmosphere
    :returns: The columns and tauross in a tuple
    """
    data = _unpack_model(fname)
    model = _loadtxt(data[23:-2])
    tauross = tauross_scale(model[:, 4], model[:, 0])
    return (model, tauross)


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

    # Some black magic we need later
    interGridShape = (2, 2, 2)
    N = reduce(mul, interGridShape)  # Multiply all values in tuple above
    idxs = [0] * N
    for cntr in range(N):
        idxs[cntr] = np.unravel_index(cntr, interGridShape)

    # Reading the models and defining the opacity intervals
    models = []
    opmin, opmax = [], []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm[0])
        opmin.append(np.amin(tatm[1]))
        opmax.append(np.amax(tatm[1]))

    # Define the grid coordinates for the interpolation
    # TODO: Put this in a function, if it is this ugly!
    try:
        c1 = [(teff-tefflow)/(teffhigh-tefflow)]
    except ZeroDivisionError:
        c1 = [0.5]
    try:
        c2 = [(logg-logglow)/(logghigh-logglow)]
    except ZeroDivisionError:
        c2 = [0.5]
    try:
        c3 = [(feh-fehlow)/(fehhigh-fehlow)]
    except ZeroDivisionError:
        c3 = [0.5]
    # Need to be like this for map_coordinates! Do not touch the line below.
    coord = np.array([c1, c2, c3])

    # Interpolate the models using the Force
    # Look at Jobovy code.
    # More optimized/clean version compared to his
    newdeck = np.zeros(models[0].shape)
    newdeck[:, 7:] = models[0][:, 7:]
    layers, columns = newdeck.shape
    for layer in range(layers):
        for column in range(columns):
            tlayer = np.zeros(interGridShape)
            for cntr in range(N):
                tlayer[idxs[cntr]] = models[cntr][layer, column]
            newdeck[layer, column] = map_coordinates(tlayer, coord)
    return newdeck


def save_model(model, params, type='kurucz95', fout='out.atm'):
    """Save the model atmosphere in the right format

    :model: The interpolated model atmosphere
    :type: Type of model atmosphere (onyl kurucz95 at the moment)
    :fout: Which place to save to
    """
    model = model[:, 0:7]
    teff, logg, feh, vt = params
    if type == 'kurucz95':
        header = 'KURUCZ\n'\
                 'Teff=%i   log g=%.2f   [Fe/H]=%.2f    vt=%.3e\n'\
                 'NTAU        %i' % (teff, logg, feh, vt, model.shape[0])
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

    _fmt = ('%15.8E', '%8.1f', '%.3E', '%.3E', '%.3E', '%.3E', '%.3E')
    np.savetxt(fout, model, header=header, footer=footer, comments='',
               delimiter=' ', fmt=_fmt)
