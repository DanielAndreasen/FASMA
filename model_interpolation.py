#!/usr/bin/python
import numpy as np
from scipy.integrate import simps
import scipy as sp
from scipy.interpolate import interp1d
import gzip

from pymoog import _get_model


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


def tauross_scale(abross, rhox, num_layers):
    """Build the tau-ross scale??

    :abross: Some model thing
    :rhox: another model thing
    :num_layers: Number of layers in the model atmosphere
    :returns: the new tau-ross scale

    """
    tauross = np.zeros(num_layers)
    # This is supposed to be the first element
    tauross[0] = abross[0] * rhox[0]
    for i in range(2, num_layers+1):
        tauross[i-1] = sp.integrate.simps(rhox[0:i], abross[0:i]) #, even='last')
        # print i, tauross[i-1]
    return tauross


def read_model(filename):
    """Read the model, return all the columns and tauross"""

    teff, num_layers = _read_header(filename)

    # This are the thermodynamical quantities.
    model = np.genfromtxt(filename, dtype=None, skiprows=23, skip_footer=2,
                          names=['RHOX', 'T', 'P', 'XNE', 'ABROSS',
                                 'ACCRAD', 'VTURB', 'tab1', 'tab2'])

    # TODO: Can this even happen? Any way, a better error would be helpful :)
    if len(model) != num_layers:
            raise Exception("FORMAT ERROR")

    # TODO: Do we need all this? Possible optimization for later
    model_rhox = model['RHOX']
    model_t = model['T']
    model_p = model['P']
    model_xne = model['XNE']
    model_abross = model['ABROSS']
    model_accrad = model['ACCRAD']
    model_vturb = model['VTURB']
    model_tab1 = model['tab1']
    model_tab2 = model['tab2']

    tauross = tauross_scale(model_abross, model_rhox, num_layers)

    return (model_rhox, model_t, model_p, model_xne, model_abross,
            model_accrad, model_vturb, tauross)




# We can also find all models in the grid. Gives back the 8 columns we want :)
models = _get_model(teff=5777, logg=4.44, feh=0.00, type='kurucz95')
tauross_all = []
model_all = []
for model in models:
    read = _unpack_model(model)
    columns = read_model(read)
    tauross = columns[-1]
    tauross_all.append(tauross)
    model_all.append(columns[0:-1])

tauross = tauross_all[0]
ntau = len(tauross)

tauross_min = min([v[-1] for v in tauross_all])
tauross_max = max([v[0] for v in tauross_all])

tauross_tmp = tauross[(tauross > tauross_max) & (tauross < tauross_min)]
f = interp1d(range(len(tauross_tmp)), tauross_tmp)
tauross_new = f(np.linspace(0, len(tauross_tmp) -1, ntau))
















# x1 = read_model('5500g45.p01')
# x2 = read_model('5750g45.p01')
# x3 = read_model('6000g45.p01')
# x4 = read_model('6250g45.p01')
# print x1[2]
# print x2[2]
# print x3[2]
# print x4[2]

# avail_models = ['5500g45.p01', '5750g45.p01', '6000g45.p01', '6250g45.p01']
# num_of_models = len(avail_models)

# # Now I will resample everything on the rosseland optical depth scale of model1
# tauross_max = []
# tauross_min = []
# # These are the min and max values of tauross
# # TODO: Here we read the model 3 times, when we only need to do it once
# # TODO: We can use numpy arrays instead of .append
# models = [read_model(model) for model in avail_models]
# for i, model in enumerate(avail_models):
#     model_i = read_model(model)
#     models.append(model_i)
#     tauross_min.append(model_i[7][0])
#     tauross_max.append(model_i[7][-1])

# max_tauross = max(tauross_max)
# min_tauross = min(tauross_min)
# print max_tauross
