#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.optimize import root, fmin_tnc, fmin_slsqp, fmin_l_bfgs_b, fmin_cobyla, minimize
from runmoog import fun_moog
from pymoog import _get_model
from model_interpolation import interpolator, save_model
import os


def callback_f(xk, convergence=None):
    os.system('test touch')
    models, nt, nl, nf = _get_model(teff=xk[0], logg=xk[1], feh=xk[2])
    n = interpolator(models, teff=(xk[0], nt), logg=(xk[1], nl), feh=(xk[2], nf))
    save_model(n, xk)


x0 = 4435, 3.54, 0.13, 2.4
bounds = [(4000, 6000), (3.00, 4.50), (-0.1, 0.1), (0.0, None)]

# a = differential_evolution(fun_moog, bounds, polish=False)
# a = fmin_l_bfgs_b(fun_moog, x0, bounds=bounds, approx_grad=True)
# a = fmin_tnc(fun_moog, x0, approx_grad=True, bounds=bounds)
# a = fmin_slsqp(fun_moog, x0 , bounds=bounds)
a = minimize(fun_moog, x0 , bounds=bounds, method='L-BFGS-B')

print a





# for i, xi in enumerate(a[0]):
    # print x0[i] - xi


# print x0
# print a[0]
# print a[1]
# print a[2]









# models, nt, nl, nf = _get_model(teff=x0[0], logg=x0[1], feh=x0[2])
# n = interpolator(models, teff=(x0[0], nt), logg=(x0[1], nl), feh=(x0[2], nf))
# save_model(n, x0)



# For a vector function
# a = root(fun_moog, x0, method='broyden1', callback=callback_f, options=opt)

