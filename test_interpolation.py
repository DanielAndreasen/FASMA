import os
import numpy as np

import pytest

from interpolation import read_model
from interpolation import interpolator
from interpolation import save_model


def test_read_model():
    fname = 'models/kurucz95/p00/5750g45.p00.gz'
    model = read_model(fname)
    assert isinstance(model, np.ndarray)
    assert model.shape[0] > model.shape[1]
    with pytest.raises(IOError):
        read_model('wrong-path')


def test_interpolator_kurucz95():
    params = (5777, 4.0, 0.04, 1.00)
    m, p = interpolator(params, save=False, atmtype='kurucz95', result=True)
    assert len(p) == len(params)
    assert list(p) == list(params)
    assert isinstance(m, np.ndarray)
    assert m.shape[0] > m.shape[1]
    with pytest.raises(ValueError):
        m, p = interpolator('params', save=False, atmtype='kurucz95', result=True)

    res = interpolator(params, save=False, atmtype='kurucz95', result=False)
    assert res is None


def test_interpolator_save():
    params = (5777, 4.0, 0.04, 1.00)
    interpolator(params)
    assert os.path.isfile('out.atm')


def test_interpolator_wrong_model():
    params = (5777, 4.0, 0.04, 1.00)
    with pytest.raises(NotImplementedError):
        interpolator(params, atmtype='wrong')


def test_save_model():
    params = (5777, 4.0, 0.04, 1.00)
    m, p = interpolator(params, save=False, atmtype='kurucz95', result=True)
    save_model(m, p, fout='test.atm')
    assert os.path.isfile('test.atm')
    os.remove('test.atm')
    with pytest.raises(NameError):
        save_model(m, p, type='wrong', fout='test.atm')
