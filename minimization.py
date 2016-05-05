#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from copy import copy


class Minimize:
    '''Minimize for best parameters given a line list with the EW method

    Inputs
    ------
    x0 : list/tuple
      The initial guess on parameters containing (Teff, logg, [Fe/H], vt)
    func : function
      The function to minimize calculating slopes for a set of parameters
    model : str
      The model atmosphere to be used
    weights : str (default: 'null')
      The weights to be used for calculating correlation slopes
    fix_teff : bool (default: False)
      Fixing the Teff when set to True
    fix_logg : bool (default: False)
      Fixing the logg when set to True
    fix_feh : bool (default: False)
      Fixing the [Fe/H] when set to True
    fix_vt : bool (default: False)
      Fixing the vt when set to True
    iterations : int (default: 160)
      Maximum number of allowed iterations
    EPcrit : float (default: 0.001)
      The criteria for convergence on abundances vs excitation potential
    RWcrit : float (default: 0.003)
      The criteria for convergence on abundances vs reduced EW
    ABdiffcrit : float (default: 0.01)
      The criteria for convergence on FeI-FeII
    MOOGv : int (default: 2014)
      The MOOG version being used
    GUI : bool (default: True)
      Change the printing behaviour according to wether it is in the GUI or not
    '''

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
        '''Get the microturbulence if this is fixed from an emperical relation
        From dwarfs the relation is from Tsantaki+ 2013, while for giants
        it is from Adibekyan+ 2015.

        Output
        ------
          Update x0[3]  --  vt
        '''
        if self.x0[1] >= 3.95:
            self.x0[3] = 6.932*self.x0[0]/10000 - 0.348*self.x0[1] - 1.437
        else:
            self.x0[3] = 3.7 - 5.1*self.x0[0]/10000


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

        Input
        -----
        i : int
          the index of the bounds we want to check. 4 values are possible
          1 means Teff
          3 means logg
          5 means [Fe/H]
          7 means vt
        '''
        if self.x0[int((i-1)/2)] < self.bounds[i-1]:
            self.x0[int((i-1)/2)] = self.bounds[i-1]
        elif self.x0[int((i-1)/2)] > self.bounds[i]:
            self.x0[int((i-1)/2)] = self.bounds[i]


    def check_convergence(self, fe_input):
        '''Check the convergence criteria for all parameters which are not fixed

        Input
        -----
        fe_input : float
          The iron abundance

        Output
        ------
        conditions : bool
          If converged return True, False otherwise
        '''
        self.slopeEP = 0.00 if self.fix_teff else self.slopeEP
        self.slopeRW = 0.00 if self.fix_vt else self.slopeRW
        self.Abdiff = 0.00 if self.fix_logg else self.Abdiff
        self.x0[2] = fe_input+7.47 if self.fix_feh else self.x0[2]

        cond1 = abs(self.slopeRW) <= self.RWcrit
        cond2 = abs(self.Abdiff) <= self.ABdiffcrit
        cond3 = abs(self.slopeEP) <= self.EPcrit
        cond4 = round(fe_input, 2) == round(self.x0[2]+7.47, 2)
        return cond1 and cond2 and cond3 and cond4


    def _bump(self, alpha):
        '''Bump to the values in the list by making a draw from a Gaussian
        distribution centered around xi with a width alpha_i

        Input
        -----
        alpha : list
          A list with the different sigma for the Gaussian distribution

        Output
        ------
          An updated parameters list
        '''
        for i, X in enumerate(zip(alpha, self.x0)):
            ai, xi = X
            sig = 0.01 if ai*xi == 0 else ai*xi
            if ai:
                self.x0[i] = np.random.normal(xi, abs(sig))


    def _format_x0(self):
        '''Format the values in x0, so first value is an integer and rest
        are floats with 2 decimals.'''
        self.x0[0] = int(self.x0[0])
        self.x0[1] = round(self.x0[1], 2)
        self.x0[2] = round(self.x0[2], 2)
        self.x0[3] = round(self.x0[3], 2)


    def minimize(self):
        '''Do the actual minimization after initializing Minimize (the class).

        Output
        ------
        x0 : list
          The final parameters
        converged : bool
          Return wether the minimization was succesful.
        '''
        step = (700, 0.50, 0.50)

        self._format_x0()
        res, self.slopeEP, self.slopeRW, abundances = self.func(self.x0, self.model, version=self.MOOGv)
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
                s = np.sign(self.slopeEP)
                step_i = s * step[0]/abs(np.log(abs(self.slopeEP)+0.0005))**3
                step_i = s*1 if abs(step_i) < 1 else step_i
                self.x0[0] += step_i
                self.check_bounds(1)

            # Step for VT
            if (abs(self.slopeRW) >= self.RWcrit) and not self.fix_vt:
                s = np.sign(self.slopeRW)
                step_i = s * step[2]/abs(np.log(abs(self.slopeRW)+0.0005))**3
                step_i = s*0.01 if abs(step_i) < 0.01 else step_i
                self.x0[3] += step_i
                self.check_bounds(7)

            # Step for logg
            if (abs(self.Abdiff) >= self.ABdiffcrit) and not self.fix_logg:
                s = -np.sign(self.Abdiff)
                step_i = s * abs(self.Abdiff)
                step_i = s*0.01 if abs(step_i) < 0.01 else step_i
                self.x0[1] += step_i
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
            res, self.slopeEP, self.slopeRW, abundances = self.func(self.x0, self.model, weights=self.weights, version=self.MOOGv)
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


def minimize_synth(p0, x_obs, y_obs, N, r, f, options):
    '''Minimize a synthetic spectrum to an observed

    Input
    -----
    p0 : list/tuple
      Initial parameters (teff, logg, feh, vt)
    x_obs : array-like
      The observed wavelengths
    y_obs : array-like
      The observed flux
    N :
    r :
    f :
    options :

    Output
    -----
    params : list
      Final parameters
    x_final : array-like
      Final wavelength
    y_final : array-like
      Final synthetic flux
    '''

    from utils import fun_moog as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec

    x = (N, r, f, x_obs, options)

    def myfunct(p, x=None, y=None, err=None):
        '''Function that return the weighted deviates
        Input
        ----
        p : list/tuple
          Parameters for the model atmosphere
        x : list
          Something something
        y : array-like
          Observed flux
        err : array-like
          Observed error

        Output
        -----
        (y-model)/err : float
          The xi-value
        '''

        #Definition of the Model spectrum to be iterated
        #N, r, f, x_obs, options = x
        x_s, y_s = func(p, atmtype=options['model'], driver='synth',
                        version=options['MOOGv'], N=N, r=r,
                        fout=f, options=options)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_s = sl(x_obs)
        model = flux_s
        status = 0
        return([status, (y-model)/err])

    #error of the flux.. this should be better calculated
    err = np.zeros(len(y_obs)) + 0.01
    fa = {'x': [N, r, f, x_obs, y_obs, options], 'y': y_obs, 'err': err}

    #set PARINFO structure
    parinfo = [{'limited': [1, 1], 'limits': [4000.0, 7500.0], 'step': 10},
               {'limited': [1, 1], 'limits': [0.5, 5.0], 'step': 0.05},
               {'limited': [1, 1], 'limits': [-2.0, 1.0], 'step': 0.01},
               {'limited': [1, 1], 'limits': [0.0, 10.0], 'step': 0.01}]

    m = mpfit(myfunct, xall=p0, functkw=fa, parinfo=parinfo)
    print('status = %s' % m.status)
    print('Iterations: %s' % m.niter)
    print('Fitted pars:%s' % m.params)
    print('Uncertainties: %s' % m.perror)

    x_s, y_s = func(m.params, atmtype=options['model'], driver='synth',
                    version=options['MOOGv'], N=N, r=r, fout=f, options=options)
    sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
    flux_final = sl(x_obs)
    #I should create a heaader with the parameters here
    save_synth_spec(x_obs, flux_final, fname='final.spec')
    return m.params, x_obs, flux_final


def mcmc_synth(x0, observed, limits):

    import emcee
    from utils import interpol_synthetic, fun_moog as func
    from utils import read_observations

    def lnlike(theta, x, y, yerr):
        teff, logg = theta
        x0 = (teff, logg, 0.0, 1.0)
        func(x0, driver='synth')
        _, _, model = interpol_synthetic(x, y, limits[0], limits[1])
        inv_sigma2 = 1.0/(yerr**2 + model**2*np.exp(2))
        return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))

    def lnprior(theta):
        teff, logg = theta
        if 5500 < teff < 6000 and 4.0 < logg < 4.9:
            return 0.0
        return -np.inf

    def lnprob(theta, x, y, yerr):
        lp = lnprior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + lnlike(theta, x, y, yerr)



    x, y = np.loadtxt(observed, unpack=True, usecols=(0, 1))
    idx = (x >= limits[0]) & (x <= limits[1])
    x, y = x[idx], y[idx]
    y /= np.median(y)
    # Normalization (use first 50 points below 1.2 as constant continuum)
    maxes = y[(y < 1.2)].argsort()[-50:][::-1]
    y /= np.median(y[maxes])


    x0 = np.array(x0)


    ndim, nwalkers = 2, 8
    Teff_step, logg_step = 50, 0.1
    pos = [x0 + np.random.randn()*np.array([Teff_step, logg_step]) for i in range(nwalkers)]

    print('starting MCMC')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, np.array([0.01]*len(x))))
    print('still doing MCMC I guess')
    sampler.run_mcmc(pos, 500)

    print('Are we done soon???')
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    print('Done!')
    # print(samples)
    print(samples.shape)
    import corner
    import matplotlib.pyplot as plt
    fig = corner.corner(samples, labels=["$Teff$", "$logg$"], truths=[5777, 4.44])
    plt.show()
