#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.optimize import minimize
from runmoog import fun_moog, fun_moog_fortran
import matplotlib.pyplot as plt


def construct_jacobian(func,epsilon):
    def jac(x, *args):
        x0 = np.asfarray(x)
        f0 = np.atleast_1d(func(*((x0,)+args)))
        jac = np.zeros([len(x0),len(f0)])
        dx = np.zeros(len(x0))
        for (i, e) in zip(range(len(x0)), epsilon):
            dx[i] = e
            jac[i] = (func(*((x0+dx,)+args)) - f0)/e
            dx[i] = 0.0
        return jac.transpose()
    return jac



Tscale = 10
x0 = 4750/Tscale, 2.54, -0.03, 1.10
bnds = ((4500/Tscale, 30000/Tscale), (1.50, 5.00), (-2.0, 0.5), (0.0, None))
cons = ({'type': 'ineq', 'fun': lambda x: 3000/Tscale < x[0] < 35000/Tscale},
        {'type': 'ineq', 'fun': lambda x: 2.5 < x[1] < 5.00},
        {'type': 'ineq', 'fun': lambda x: -3.0 < x[2] < 1.0},
        {'type': 'ineq', 'fun': lambda x: 0.0 < x[3] < 6.00})

## For COBYLA
opt = {'tol': 1e-8,
       # 'rhobeg': 1.0,
        }
# a = minimize(fun_moog_fortran, x0, method='COBYLA', options=opt)
a = minimize(fun_moog, x0, method='COBYLA', options=opt)







# a = minimize(fun_moog_fortran, x0,
#              jac=construct_jacobian(fun_moog_fortran, [5, 1e-2, 1e-2, 1e-2]),
#              bounds=bnds, method='SLSQP')

# a = minimize(fun_moog_fortran, x0, method='Nelder-Mead')


## For TNC
# opt = {'ftol': 0.05,
#        'gtol': 0.01,
#        'scale': [0.001, 1.0, 1.0, 1.0],
#        # 'rescale': 0,
#        'eps': 0.01
#         }



# a = minimize(fun_moog, x0 , bounds=bnds, method='TNC', options=opt)

# a = differential_evolution(fun_moog, bnds, polish=False)
# a = minimize(fun_moog, x0 , bounds=bnds, method='L-BFGS-B', options=opt)

# Does not work
# a = minimize(fun_moog, x0 , bounds=bnds, method='SLSQP', options=opt)

print a
print ''
print 'Temperature: %i' % (float(a.x[0])*Tscale)
print '      log g: %.2f' % a.x[1]
print '     [Fe/H]: %.2f' % a.x[2]
print '         vt: %.2f' % a.x[3]
