#!/usr/bin/python
from __future__ import division
import numpy as np
from operator import mul
from scipy.ndimage import _nd_image
import gzip
from scipy import ndimage
from scipy import interpolate


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
    data = np.zeros((row, col))
    for rowi, line in enumerate(lines):
        try:
            data[rowi, :] = map(float, line.split())
        except ValueError:
            return data[0:rowi, :]
    return data


def read_model(fname):
    """Read the model atmosphere

    :fname: The gz file of the model atmosphere
    :returns: The columns and tauross in a tuple
    """
    data = _unpack_model(fname)
    model = _loadtxt(data[23:-1])
    tauross = tauross_scale(model[:, 4], model[:, 0])
    return (model, tauross)


def interpolator(mnames, teff, logg, feh):
    """The function to call from the main program (MOOGme.py)

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
    c1 = [(teff-tefflow)/(teffhigh-tefflow)]
    c2 = [(logg-logglow)/(logghigh-logglow)]
    c3 = [(feh-fehlow)/(fehhigh-fehlow)]

    # Need to be like this for map_coordinates! Do not touch the line below.
    coord = np.array([c1, c2, c3])

    # Interpolate the models using the Force
    # Look at Jobovy code.
    # More optimized/clean version compared to his

    layers = min([model.shape[0] for model in models])
    columns = min([model.shape[1] for model in models])
    newdeck = np.zeros((layers, columns))
    # newdeck[:, 7:] = models[0][:, 7:]
    for layer in range(layers):
        for column in range(columns):
            tlayer = np.zeros(interGridShape)
            for cntr in range(N):
                tlayer[idxs[cntr]] = models[cntr][layer, column]
            # newdeck[layer, column] = map_coordinates(tlayer, coord)
            newdeck[layer, column] = ndimage.interpolation.map_coordinates(tlayer, coord, order=1)
    return newdeck

def interpolatorN7(mnames, teff, logg, feh):
    """The function to call from the main program (MOOGme.py)

    :mnames: As generated from _get_models
    :teff: Requested Effective Temperature and the four closest models in
           the gridpoints
    :logg: Requested Surface gravity and the two closest models
    :feh: Requested metallicity and the two closest models
    :out: The interpolated model saved in this file
    """

    def distance_vector(teff, logg, feh, teffvector, loggvector, fehvector):
        """This function computes the distance_vector
        :teff: Requested Effective Temperature
        :logg: Requested logg
        :feh:  requested feh
        :teffvector: vector with closest 4 teff
        :loggvector: vector with closest 2 logg
        :fehvector: vector with closest 2 Fe/H
        :vector: output vector with the distance of each gridpoint to the
        requested parameters.log10(
        """
        vector = []
        for element in teffvector:
            for element2 in loggvector:
                for element3 in fehvector:
                    # x1 = np.log10(teff)-np.log10(element)
                    x1 = (teff-element)/teff
                    x2 = logg-element2
                    x3 = feh-element3
                    vector.append(x1 + x2/logg + x3)
        return vector


    def dist(teff, logg, feh, mname):
        m = mname.rpartition('/')[-1]
        m = m.split('.')
        tm = int(m[0][:-3])
        gm = float(m[0][-2:])/10
        if m[1].startswith('m'):
            fm = -float(m[1][1:])/10
        else:
            fm = float(m[1][1:])/10
        x1 = (teff-tm)/teff
        x2 = logg-gm
        x3 = feh-fm
        return 30*x1 + 1.5*x2 + x3


    # TODO: This is super ugly way of doing stuff!
    # teff, tefflow1, tefflow2, teffhigh1, teffhigh2 = teff[0], teff[1][0], teff[1][1], teff[1][2], teff[1][3]
    teff, tefflow1, teffhigh1 = teff[0], teff[1][0], teff[1][1]
    logg, logglow, logghigh = logg[0], logg[1][0], logg[1][1]
    feh, fehlow, fehhigh = feh[0], feh[1][0], feh[1][1]

    # teff_vector = [tefflow1, tefflow2, teffhigh1, teffhigh2]
    teff_vector = [tefflow1, teffhigh1]
    logg_vector = [logglow, logghigh]
    feh_vector  = [fehlow, fehhigh]

    distance = distance_vector(teff, logg, feh, teff_vector, logg_vector, feh_vector)

    from matplotlib import pyplot as plt

    # Reading the models and defining the opacity intervals
    models = []
    distance = []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm[0])
        distance.append(dist(teff, logg, feh, mname))
    idx = np.argsort(distance)
    distance = list(np.array(distance)[idx])
    models = list(np.array(models)[idx])

    layers = min([model.shape[0] for model in models])
    columns = min([model.shape[1] for model in models])
    newatm = np.zeros((layers, columns))
    # newdeck[:, 7:] = models[0][:, 7:]
    for layer in range(layers):
        for column in range(columns):
            parameter_vector = []
            for model in models:
                parameter_vector.append(model[layer,column])
            f = interpolate.interp1d(distance, parameter_vector, kind='cubic')
            # plt.plot(distance, parameter_vector, 'o-k')
            # x = np.linspace(distance[0], distance[-1], 100)
            # plt.plot(x, f(x), '-r')
            # plt.plot([0], f(0), '-or')
            # plt.show()
            newatm[layer,column] = f(0)
    return newatm#, models[0]







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
                 'Teff= %i   log g= %.2f\n'\
                 'NTAU        %i' % (teff, logg, model.shape[0])
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
             '         8.1    822.0     22.1' % (vt*1e5, feh, 7.47+feh)

    _fmt = ('%15.8E', '%8.1f', '%.3E', '%.3E', '%.3E', '%.3E', '%.3E')
    np.savetxt(fout, model, header=header, footer=footer, comments='',
               delimiter=' ', fmt=_fmt)



if __name__ == '__main__':
    from utils import _get_model
    import matplotlib.pyplot as plt
    teff = 5250
    logg = 4.2
    feh = 0.0
    mnames, teffmod, loggmod, fehmod = _get_model(teff, logg, feh, atmtype='kurucz95')
    new_atm, model = interpolatorN7(mnames, [teff, teffmod], [logg,loggmod], [feh,fehmod])
    # for i in range(1, 9):
    #     plt.plot(new_atm[:,0], (new_atm[:,i]-model[:,i])/model[:,i] * 100, 'ro')
    #     plt.show()
