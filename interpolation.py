#!/usr/bin/python
from __future__ import division
import numpy as np
import gzip
from scipy.interpolate import griddata


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


def read_model(fname):
    """Read the model atmosphere

    :fname: The gz file of the model atmosphere
    :returns: The columns and tauross in a tuple
    """
    f = gzip.open(fname, compresslevel=1)
    data = f.readlines()
    model = np.loadtxt(data[23:-2])
    return model


def interpolator(mnames, teff, logg, feh):
    """This is a new approach based on a scipy interpolator.
    Resembles the original interpolator we used but with a change

    :mnames: As generated from GetModels
    :teff: Requested Effective Temperature and the two closest models in
           the gridpoints
    :logg: Requested Surface gravity and the two closest models
    :feh: Requested metallicity and the two closest models
    :out: The interpolated model saved in this file
    """

    # Making the surrounding grid points
    gridpoints = []
    for temp in teff[1]:
        for grav in logg[1]:
            for metal in feh[1]:
                gridpoints.append((temp,grav,metal))
    gridpoints = np.asarray(gridpoints)

    # Define the points to obtain at the end
    teff = teff[0]
    logg = logg[0]
    feh = feh[0]

    # Reading the models
    models = []
    for mname in mnames:
        tatm = read_model(mname)
        models.append(tatm)

    layers = min([model.shape[0] for model in models])
    columns = min([model.shape[1] for model in models])
    newatm = np.zeros((layers, columns))
    for layer in range(layers):
        for column in range(columns):
            tlayer = np.zeros(len(models))
            for idx, model in enumerate(models):
                tlayer[idx] = model[layer, column]
            newatm[layer, column] = griddata(gridpoints, tlayer, (teff, logg, feh), method='linear', rescale=True)
    return newatm


def save_model(model, params, type='kurucz95', fout='out.atm'):
    """Save the model atmosphere in the right format

    :model: The interpolated model atmosphere
    :params: Atmospheric parameters in usual order
    :type: Type of model atmosphere
    :fout: Which place to save to
    """
    model = model[:, 0:7]
    teff, logg, feh, vt = params
    if type == 'kurucz95' or type == 'apogee_kurucz' or type == 'marcs':
        header = 'KURUCZ\n'\
                 'Teff= %i   log g= %.2f\n'\
                 'NTAU        %i' % (teff, logg, model.shape[0])
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
    np.savetxt(fout, model, header=header, footer=footer, comments='', delimiter=' ', fmt=_fmt)
