#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from model_interpolation import interpolator, save_model
from utils import _get_model

teff = 6454
logg = 4.0

models, nt, nl, nf = _get_model(teff=teff, logg=logg, feh=0.04)
inter_model = interpolator(models, teff=(teff, nt), logg=(logg, nl), feh=(0.04, nf))
save_model(inter_model, params=(teff, logg, 0.04, 1.0))
