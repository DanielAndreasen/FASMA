#!/usr/bin/python
from __future__ import division
import numpy as np
import gzip
from scipy.interpolate import griddata
from utils import GetModels


def read_model(fname):
    '''Read the model atmosphere

    Input
    -----
    fname : str
      The gz file of the model atmosphere.

    Output
    ------
    model : ndarray
      The correct atmosphere, the columns and tauross in a tuple
    '''
    f = gzip.open(fname, compresslevel=1)
    data = f.readlines()
    model = np.loadtxt(data[23:-2])
    return model


def interpolator(params, save=True, atmtype='kurucz95'):
    '''This is a new approach based on a scipy interpolator.
    Resembles the original interpolator we used but with a change

    Input
    -----
    params : list of length 3
      Teff, logg, [Fe/H] desired.
    save : bool
      Wether the new atmosphere should be saved. Default is True.
    atmtype : str
      The atmosphere models being used. Default is Kurucz95.

    Output
    ------
    newatm : ndarray
      New interpolated atmosphere.
    '''

    m = GetModels(params[0], params[1], params[2], atmtype=atmtype)
    mdict = m.getmodels()
    mnames = mdict['models']
    teff = mdict['teff']
    logg = mdict['logg']
    feh = mdict['feh']

    # Making the surrounding grid points
    gridpoints = []
    for temp in teff[1]:
        for grav in logg[1]:
            for metal in feh[1]:
                gridpoints.append((temp, grav, metal))
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
    if save:
        save_model(newatm, params, type=atmtype)
    return newatm


def save_model(model, params, type='kurucz95', fout='out.atm'):
    '''Save the model atmosphere in the right format

    Input
    -----
    model : ndarray
      The interpolated model atmosphere.
    params : list
      Teff, logg, [Fe/H], vt of the interpolated atmosphere.
    type : str
      Type of atmospheric parameters. Default is Kurucz95
    fout : str
      Name of the saved atmosphere. Default is out.atm

    Output
    ------
    Atmospheric model.
    '''
    model = model[:, 0:7]
    teff, logg, feh, vt = params
    if type in ['kurucz95', 'apogee_kurucz', 'marcs']:
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
