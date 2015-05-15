#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from scipy.optimize import minimize
from runmoog import fun_moog

Tscale = 1000
x0 = 4075/Tscale, 2.04, -0.23, 0.4
bnds = ((4500/Tscale, 30000/Tscale), (1.50, 5.00), (-2.0, 0.5), (0.0, None))
cons = ({'type': 'ineq', 'fun': lambda x: 3000/Tscale < x[0] < 35000/Tscale},
        {'type': 'ineq', 'fun': lambda x: 2.5 < x[1] < 5.00},
        {'type': 'ineq', 'fun': lambda x: -3.0 < x[2] < 1.0},
        {'type': 'ineq', 'fun': lambda x: 0.0 < x[3] < 6.00})

## For TNC
# opt = {'ftol': 0.05,
#        'gtol': 0.01,
#        'scale': [0.001, 1.0, 1.0, 1.0],
#        # 'rescale': 0,
#        'eps': 0.01
#         }


## For COBYLA
opt = {'tol': 0.05,
        'rhobeg': [2.1, 1.0, 0.1, 0.1],
        }

# a = minimize(fun_moog, x0 , bounds=bnds, method='TNC', options=opt)
a = minimize(fun_moog, x0, method='COBYLA', constraints=cons, options=opt)

# a = differential_evolution(fun_moog, bnds, polish=False)
# a = minimize(fun_moog, x0 , bounds=bnds, method='L-BFGS-B', options=opt)

# Does not work
# a = minimize(fun_moog, x0 , bounds=bnds, method='Nelder-Mead', options=opt)
# a = minimize(fun_moog, x0 , bounds=bnds, method='SLSQP', options=opt)

print a
print ''
print 'Temperature: %i' % (float(a.x[0])*Tscale)
print '      log g: %.2f' % a.x[1]
print '     [Fe/H]: %.2f' % a.x[2]
print '         vt: %.2f' % a.x[3]
