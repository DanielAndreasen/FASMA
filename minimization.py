#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from copy import copy


class Minimize:
    '''Minimize for best parameters given a line list'''

    def __init__(self, x0, func, model, weights='null',
                 fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
                 iterations=160, EPcrit=0.001, RWcrit=0.003, ABdiffcrit=0.01,
                 MOOGv=2014, GUI=True, **kwargs):
        self.x0 = x0
        self.func = func
        self.model = model
        self.weights = weights
        self.fix_teff = fix_teff
        self.fix_logg = fix_logg
        self.fix_feh = fix_feh
        self.fix_vt = fix_vt
        self.maxiterations = iterations
        self.iteration = 0
        self.EPcrit = EPcrit
        self.RWcrit = RWcrit
        self.ABdiffcrit = ABdiffcrit
        self.MOOGv = MOOGv
        self.GUI = GUI
        if self.model.lower() == 'kurucz95':
            self.bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]
        if self.model.lower() == 'apogee_kurucz':
            self.bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99]
        if self.model.lower() == 'marcs':
            self.bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99]

    def _getMic(self):
        '''Get the microturbulence if this is fixed'''
        if self.x0[1] >= 3.95:
            self.x0[3] = 6.932*self.x0[0]/10000 - 0.348*self.x0[1] - 1.437
        else:
            self.x0[3] = 2.72 - 0.457*self.x0[1] + 0.072*self.x0[2]

    def print_format(self):
        '''Print the stellar atmospheric parameters in a nice format'''
        rest = self.x0 + list((self.slopeEP, self.slopeRW, self.Abdiff))
        if self.iteration == 0:
            if self.GUI:
                print(' i     Teff       logg     [Fe/H]    vt    EPslope    RWslope    |FeI-FeII|')
                print('-' * 99)
            else:
                print(' i    Teff    logg    [Fe/H]    vt    EPslope    RWslope    |FeI-FeII|')
                print('-' * 70)
        else:
            print '{:4d}{:>6d}{:>8.2f}{:>+9.2f}{:>8.2f}{:>+9.3f}{:>+11.3f}{:>11.2f}'.format(self.iteration, *rest)

    def check_bounds(self, i):
        '''
        Function which checks if parameters are within bounds.
        Input - parameter: what we want to check; bounds: ze bounds;
        i: the index of the bounds we want to check'''
        if self.x0[int((i-1)/2)] < self.bounds[i-1]:
            self.x0[int((i-1)/2)] = self.bounds[i-1]
        elif self.x0[int((i-1)/2)] > self.bounds[i]:
            self.x0[int((i-1)/2)] = self.bounds[i]

    def check_convergence(self, fe_input):
        '''Check the convergence criteria'''
        self.slopeEP = 0.00 if self.fix_teff else self.slopeEP
        self.slopeRW = 0.00 if self.fix_vt else self.slopeRW
        self.Abdiff = 0.00 if self.fix_logg else self.Abdiff
        fe_input = self.x0[2]+7.47 if self.fix_feh else fe_input

        cond1 = abs(self.slopeRW) <= self.RWcrit
        cond2 = abs(self.Abdiff) <= self.ABdiffcrit
        cond3 = abs(self.slopeEP) <= self.EPcrit
        cond4 = round(fe_input, 2) == round(self.x0[2]+7.47, 2)
        return cond1 and cond2 and cond3 and cond4

    def _bump(self, alpha):
        '''Bump to the values in the list, x'''
        for i, X in enumerate(zip(alpha, self.x0)):
            ai, xi = X
            sig = 0.01 if ai*xi == 0 else ai*xi
            if ai:
                self.x0[i] = np.random.normal(xi, abs(sig))

    def _format_x0(self):
        '''Format the values in x0, so first value is an integer'''
        self.x0[0] = int(self.x0[0])
        self.x0[1] = round(self.x0[1], 2)
        self.x0[2] = round(self.x0[2], 2)
        self.x0[3] = round(self.x0[3], 2)

    def minimize(self):
        self._format_x0()
        res, self.slopeEP, self.slopeRW, abundances, self.x0 = self.func(self.x0, self.model, version=self.MOOGv)
        self.Abdiff = np.diff(abundances)[0]
        self.x0 = list(self.x0)

        if self.check_convergence(abundances[0]):
            return self.x0, True

        parameters = [copy(self.x0)]
        best = {}
        # Print the header before starting
        self.print_format()

        while self.iteration < self.maxiterations:
            # Step for Teff
            if (abs(self.slopeEP) >= self.EPcrit) and not self.fix_teff:
                self.x0[0] += 2000*self.slopeEP
                self.check_bounds(1)

            # Step for VT
            if (abs(self.slopeRW) >= self.RWcrit) and not self.fix_vt:
                self.x0[3] += 1.5*self.slopeRW
                self.check_bounds(7)

            # Step for logg
            if (abs(self.Abdiff) >= self.ABdiffcrit) and not self.fix_logg:
                self.x0[1] -= self.Abdiff
                self.check_bounds(3)

            # Step for [Fe/H]
            if not self.fix_feh:
                self.x0[2] = abundances[0]-7.47
                self.check_bounds(5)

            if self.fix_vt:
                self._getMic()  # Reset the microturbulence
                self.check_bounds(7)

            if self.x0 in parameters:
                alpha = [0] * 4
                alpha[0] = abs(self.slopeEP) if not self.fix_teff else 0
                alpha[1] = abs(self.Abdiff) if not self.fix_logg else 0
                alpha[2] = 0.01 if not self.fix_feh else 0
                alpha[3] = abs(self.slopeRW) if not self.fix_vt else 0
                self._bump(alpha)
                self.check_bounds(1)
                self.check_bounds(3)
                self.check_bounds(5)
                self.check_bounds(7)
            parameters.append(copy(self.x0))

            self._format_x0()
            res, self.slopeEP, self.slopeRW, abundances, self.x0 = self.func(self.x0, self.model, weights=self.weights, version=self.MOOGv)
            self.Abdiff = np.diff(abundances)[0]
            self.iteration += 1
            self.print_format()
            best[res] = parameters[-1]
            if self.check_convergence(abundances[0]):
                print '\nStopped in %i iterations' % self.iteration
                return self.x0, True

        print '\nStopped in %i iterations' % self.iteration
        if self.check_convergence(abundances[0]):
            return self.x0, True
        else:
            # Return the best solution rather than the last iteration
            _ = self.func(best[min(best.keys())], self.model, weights=self.weights, version=self.MOOGv)
            return best[min(best.keys())], False
