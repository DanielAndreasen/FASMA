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


    def _interpol(models, param, kind, layers, columns):
        """
        models:  dictionary of the models
        param: The value we try to interpolate to
        kind: The parameter: 'teff', 'logg', 'feh'
        layers: The number of layers in the models
        columns: The number of columns in the models
        """
        newatm = np.zeros((layers, columns))
        N = len(models.keys())
        for layer in range(layers):
            for column in range(columns):
                dist = np.zeros(N)
                params = np.zeros(N)
                i = 0
                for k, v in models.iteritems():
                    if kind == 'teff':
                        kv = int(k.split('g')[0])
                    elif kind == 'logg':
                        kv = float(k[:2])/10
                    else:
                        if k[-3] == 'p':
                            kv = float(k[-2:])/10
                        else:
                            kv = -float(k[-2:])/10
                    dist[i] = kv-param
                    params[i] = v[layer, column]
                    i += 1
                # Sort the distance and parameter vector
                idx = np.argsort(dist)
                dist = dist[idx]
                params = params[idx]
                if N > 2:  # "interpolation" of temperature
                    p = np.polyfit(dist, params, deg=2)
                    p0 = np.polyval(p, 0)
                else:  # "interpolation" of logg and [Fe/H]
                    p = np.polyfit(dist, params, deg=1)
                    p0 = np.polyval(p, 0)
                newatm[layer,column] = p0
        if kind == 'teff':
            name = k.strip('.gz').split('g')[-1]
        else:
            name = k
        return {name: newatm}



    # TODO: This is super ugly way of doing stuff!
    teff, teff1, teff2, teff3, teff4 = teff[0], teff[1][0], teff[1][1], teff[1][2], teff[1][3]
    logg, logg1, logg2 = logg[0], logg[1][0], logg[1][1]
    feh, feh1, feh2 = feh[0], feh[1][0], feh[1][1]
    teff_vector = [teff1, teff2, teff3, teff4]
    logg_vector = [logg1, logg2]
    feh_vector  = [feh1, feh2]

    models = dict()
    for mname in mnames:
        tatm = read_model(mname)
        models[mname.rpartition('/')[2]] = tatm

    layers = min([model.shape[0] for model in models.values()])
    columns = min([model.shape[1] for model in models.values()])

    # Interpolate over Teff first
    models_t1 = dict()
    models_t2 = dict()
    models_t3 = dict()
    models_t4 = dict()
    for model in models.keys():
        if 'g'+str(logg1).replace('.', '') in model and str(feh1).replace('.', '').replace('-','')+'.' in model:
            models_t1[model] = models[model]
        elif 'g'+str(logg1).replace('.', '') in model and str(feh2).replace('.', '').replace('-','')+'.' in model:
            models_t2[model] = models[model]
        elif 'g'+str(logg2).replace('.', '') in model and str(feh1).replace('.', '').replace('-','')+'.' in model:
            models_t3[model] = models[model]
        elif 'g'+str(logg2).replace('.', '') in model and str(feh2).replace('.', '').replace('-','')+'.' in model:
            models_t4[model] = models[model]

    t1 = _interpol(models_t1, teff, 'teff', layers, columns)
    t2 = _interpol(models_t2, teff, 'teff', layers, columns)
    t3 = _interpol(models_t3, teff, 'teff', layers, columns)
    t4 = _interpol(models_t4, teff, 'teff', layers, columns)


    # Interpolate Fe/H
    L1 = t1
    L2 = dict()

    if L1.keys()[0][-3:] == t2.keys()[0][-3:]:
        L1[t2.keys()[0]] = t2.values()[0]
    else:
        L2[t2.keys()[0]] = t2.values()[0]

    if L1.keys()[0][-3:] == t3.keys()[0][-3:]:
        L1[t3.keys()[0]] = t3.values()[0]
    else:
        L2[t3.keys()[0]] = t3.values()[0]

    if L1.keys()[0][-3:] == t4.keys()[0][-3:]:
        L1[t4.keys()[0]] = t4.values()[0]
    else:
        L2[t4.keys()[0]] = t4.values()[0]

    l1 = _interpol(L1, logg, 'logg', layers, columns)
    l2 = _interpol(L2, logg, 'logg', layers, columns)

    # Interpolate logg
    M = l1
    M.update(l2)
    newatm = _interpol(M, feh, 'feh', layers, columns)

    return newatm.values()[0]


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


def interpolator_scipy(mnames, teff, logg, feh):
    """This is a new approach based on a scipy interpolator. Resembles the
    original interpolator we used but with a change

    :mnames: As generated from _get_models
    :teff: Requested Effective Temperature and the two closest models in
           the gridpoints
    :logg: Requested Surface gravity and the two closest models
    :feh: Requested metallicity and the two closest models
    :out: The interpolated model saved in this file
    """

    #MAKING THE GRIDPOINTS
    gridpoints =[]
    for temp in teff[1]:
        for grav in logg[1]:
            for metal in feh[1]:
                gridpoints.append((temp,grav,metal))
    gridpoints = np.asarray(gridpoints)

    #DEFINE THE POINTS TO OBTAIN IN THE END
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
    newdeck = np.zeros((layers, columns))
    for layer in range(layers):
        for column in range(columns):
            tlayer = np.zeros(len(models))
            for cntr in range(len(models)):
                tlayer[cntr] = models[cntr][layer, column]
            newdeck[layer, column] = griddata(gridpoints, tlayer, (teff, logg, feh), method='linear', rescale=True)
    return newdeck


if __name__ == '__main__':
    from utils import _get_model
    teff = 5300
    logg = 4.2
    feh = 0.02
    mnames, teffmod, loggmod, fehmod = _get_model(teff, logg, feh, atmtype='kurucz95')
    new_atm = interpolator_scipy(mnames, [teff, teffmod], [logg,loggmod], [feh,fehmod])
    # new_atml = interpolatorN7(mnames, [teff, teffmod], [logg,loggmod], [feh,fehmod])
