import os
import numpy as np
import pandas as pd

import pytest

from utils import kurucz95
from utils import GetModels
from utils import fun_moog
from utils import Readmoog
from utils import error
from utils import slope
from utils import _update_par

np.random.seed(42)


def test_getModels():
    teff, logg, feh = 5777, 4.44, 0.00
    atmtype = 'kurucz95'
    # Init stuff
    m = GetModels(teff, logg, feh, atmtype)
    assert m.teff == teff
    assert m.logg == logg
    assert m.feh == feh
    assert m.atmtype == atmtype
    for atmtype in ('kurucz95', 'apogee_kurucz', 'marcs'):
        with pytest.raises(ValueError):
            m = GetModels(50, logg, feh, atmtype)
        with pytest.raises(ValueError):
            m = GetModels(500000, logg, feh, atmtype)
        with pytest.raises(ValueError):
            m = GetModels(teff, -10, feh, atmtype)
        with pytest.raises(ValueError):
            m = GetModels(teff, 10, feh, atmtype)
        with pytest.raises(ValueError):
            m = GetModels(teff, logg, -10, atmtype)
        with pytest.raises(ValueError):
            m = GetModels(teff, logg, 10, atmtype)

    # neighbour method
    n = m.neighbour(kurucz95['teff'], teff, k=4)
    assert len(n) == 4
    assert n[1] <= teff <= n[2]
    n = m.neighbour(kurucz95['logg'], logg, k=2)
    assert len(n) == 2
    assert n[0] <= logg <= n[1]
    n = m.neighbour(kurucz95['feh'], feh, k=2)
    assert len(n) == 2
    assert n[0] <= feh <= n[1]

    # getmodels method
    m = m.getmodels()
    assert isinstance(m, dict)
    assert isinstance(m['teff'], tuple)
    assert isinstance(m['teff'][0], int)
    assert isinstance(m['teff'][1], list)
    assert isinstance(m['teff'][1][0], int)
    assert isinstance(m['logg'], tuple)
    assert isinstance(m['logg'][0], float)
    assert isinstance(m['logg'][1], list)
    assert isinstance(m['logg'][1][0], float)
    assert isinstance(m['feh'], tuple)
    assert isinstance(m['feh'][0], float)
    assert isinstance(m['feh'][1], list)
    assert isinstance(m['feh'][1][0], float)
    assert isinstance(m['models'], list)
    assert isinstance(m['models'][0], str)

    # Find the gaps in the grid
    m = GetModels(8001, 1.01, -1.01, 'kurucz95').getmodels()
    assert isinstance(m, dict)

    # Get to the edge of the grid
    teff = 39000
    m = GetModels(teff, logg, feh, 'kurucz95')
    m = m.getmodels()
    assert len(m['teff'][1]) == 2

    # Get 'error' on logg
    m = GetModels(8001, 1.01, 0.0, 'kurucz95').getmodels()
    assert isinstance(m, dict)


def test_fun_moog():
    p = (5777, 4.44, 0.00, 1.00)
    res, EPs, RWs, abundances, x = fun_moog(p, 'kurucz95')
    assert isinstance(res, float)
    assert isinstance(EPs, float)
    assert isinstance(RWs, float)
    assert len(abundances) == 2
    assert isinstance(abundances, list)
    assert isinstance(abundances[0], float)
    assert isinstance(abundances[1], float)
    assert len(x) == len(p)
    assert isinstance(x, list)
    assert isinstance(x[0], int)
    assert isinstance(x[1], float)
    assert isinstance(x[2], float)
    assert isinstance(x[3], float)
    assert list(p) == list(x)
    assert os.path.isfile('summary.out')


def test_Readmoog():
    # Be sure to have some data to analyse
    p = (5777, 4.44, 0.00, 1.00)
    # res, EPs, RWs, abundances, x = fun_moog(p, 'kurucz95')

    # Init stuff
    r = Readmoog()
    assert list(p) == list(r.params)
    r = Readmoog(params=p)
    assert r.parameters() == p

    # Fe statistic
    f = r.fe_statistics()
    assert len(f) == 8
    assert isinstance(f[0], float)
    assert isinstance(f[1], float)
    assert isinstance(f[2], float)
    assert isinstance(f[3], float)
    assert isinstance(f[4], float)
    assert isinstance(f[5], float)
    assert isinstance(f[6], np.ndarray)
    assert isinstance(f[7], np.ndarray)
    assert f[6].shape[0] > f[7].shape[0]
    assert f[6].shape[1] == f[7].shape[1]

    # elements method
    e = r.elements()
    assert isinstance(e[0], list)
    assert isinstance(e[0][0], str)
    assert isinstance(e[1], list)
    assert isinstance(e[1][0], float)

    # all_table method
    t = r.all_table()
    assert isinstance(t, pd.DataFrame)
    assert t.shape[0] > t.shape[1]
    assert len(t.columns) <= 9
    cols = 'wavelength,EP,logGF,EWin,logRWin,abund,delavg'.split(',')
    for col in cols:
        assert col in t.columns

    # atomNameFromMOOG method
    assert r.atomNameFromMOOG('26.0') == 'FeI'
    assert r.atomNameFromMOOG('26.1') == 'FeII'
    assert r.atomNameFromMOOG('6.2') == 'CIII'
    assert r.atomNameFromMOOG('42.0') == 'MoI'


def test_error():
    ll = 'sun_harps_ganymede.moog'
    p0 = (5777, 4.44, 0.00, 1.00)
    p = error(ll, True, p0, 'kurucz95')
    for i in range(4):
        assert p0[i] == p[2*i]
    assert isinstance(p[1], int)
    assert isinstance(p[3], float)
    assert isinstance(p[5], float)
    assert isinstance(p[7], float)
    p = error(ll, False, p0, 'kurucz95')
    for i in range(4):
        assert p0[i] == p[2*i]
    assert isinstance(p[1], int)
    assert isinstance(p[3], float)
    assert isinstance(p[5], float)
    assert isinstance(p[7], float)


def test_slope():
    data = [[i for i in range(10)],
            [2*i+np.random.rand()*0.1 for i in range(10)]]
    for weight in ('null', 'sigma', 'mad', 'wrong'):
        a, b = slope(data, weights=weight)
        assert a - 2.0 < 0.1
        assert min(b) > 0
        assert max(b) - 1 < 0.1


def test_update_par():
    with pytest.raises(IOError):
        _update_par(line_list='wrong-file.moog')

    ll = 'linelist/sun_harps_ganymede.moog'
    _update_par(line_list=ll)
    assert os.path.isfile('batch.par')
    _update_par(line_list=ll, lines=0)
    assert os.path.isfile('batch.par')
    _update_par(line_list=ll, plotpars=True)
    assert os.path.isfile('batch.par')
