#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
from copy import copy
import time


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
        self.x0[2] = fe_input+7.47 if not self.fix_feh else self.x0[2]

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


class Minimize_synth:
    '''Minimize a synthetic spectrum to an observed

    Input
    -----
    p0 : list
      Initial parameters (teff, logg, feh, vt)
    x_obs : ndarray
      Observed wavelength
    y_obs : ndarray
      Observed flux
    r : ndarray
      ranges of the intervals
    fout : str
      Input line list files

    Output
    -----
    params : list
      Final parameters
    x_final : ndarray
      Final wavelength
    y_final : ndarray
      Final synthetic flux
    '''

    from utils import fun_moog_synth as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec

    def __init__(self, p0, x_obs, y_obs, r, fout, model='kurucz95',
                 fix_teff=None, fix_logg=None, fix_feh=None, fix_vt=None,
                 fix_vmac=None, fix_vsini=None, **kwargs):

        self.p0 = p0
        self.x_obs = x_obs
        self.y_obs = y_obs
        self.r = r
        self.fout = fout
        self.model = model
        self.fix_teff = 1 if fix_teff else 0
        self.fix_logg = 1 if fix_logg else 0
        self.fix_feh = 1 if fix_feh else 0
        self.fix_vt = 1 if fix_vt else 0
        self.fix_vmac = 1 if fix_vmac else 0
        self.fix_vsini = 1 if fix_vsini else 0

        # Setting up the bounds
        if self.model.lower() == 'kurucz95':
            self.bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99, 0, 50, 0, 100]
        if self.model.lower() == 'apogee_kurucz':
            self.bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99, 0, 50, 0, 100]
        if self.model.lower() == 'marcs':
            self.bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99, 0, 50, 0, 100]

        # Setting up PARINFO for mpfit
        teff_info  = {'limited': [1, 1], 'limits': self.bounds[0:2],   'step': 30,   'mpside': 2, 'fixed': fix_teff}
        feh_info   = {'limited': [1, 1], 'limits': self.bounds[4:6],   'step': 0.05, 'mpside': 2, 'fixed': fix_feh}
        logg_info  = {'limited': [1, 1], 'limits': self.bounds[2:4],   'step': 0.2,  'mpside': 2, 'fixed': fix_logg}
        vt_info    = {'limited': [1, 1], 'limits': self.bounds[6:8],   'step': 0.3,  'mpside': 2, 'fixed': fix_vt}
        vmac_info  = {'limited': [1, 1], 'limits': self.bounds[8:10],  'step': 0.5,  'mpside': 2, 'fixed': fix_vmac}
        vsini_info = {'limited': [1, 1], 'limits': self.bounds[10:12], 'step': 0.5,  'mpside': 2, 'fixed': fix_vsini}
        self.parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

        # Setting up keyword arguments for myfunct
        self.fa = {'x_obs': x_obs, 'r': r, 'fout': fout, 'model': model,
                   'y': y_obs, 'options': kwargs}

    def bounds(self, i, p):
        if p[int((i-1)/2)] < self.bounds[i-1]:
            p[int((i-1)/2)] = self.bounds[i-1]
        elif p[int((i-1)/2)] > self.bounds[i]:
            p[int((i-1)/2)] = self.bounds[i]
        return p

    def myfunct(self, p, y=None, **kwargs):
        '''Function that return the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters for the model atmosphere
        x_obs : ndarray
          Wavelength
        r : ndarray
          ranges of the intervals
        fout : str
          Line list file
        model : str
          Model atmosphere type
        y : ndarray
          Observed flux


        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observation
        '''

        options = kwargs['options']
        for i in range(1, 12, 2):
            p = self.bounds(i, p)

        x_s, y_s = func(p, atmtype=model, driver='synth', r=self.r, fout=self.fout, **options)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_s = sl(self.x_obs)
        ymodel = flux_s
        # Error on the flux #needs corrections
        err = np.zeros(len(y)) + 0.01
        status = 0
        return([status, (y-ymodel)/err])

    def minimize(self):
        start_time = time.time()
        m = mpfit(self.myfunct, xall=self.p0, parinfo=self.parinfo,
                  ftol=1e-5, xtol=1e-5, gtol=1e-10, functkw=self.fa)
        end_time = time.time()-start_time
        print('status = %s' % m.status)
        print('Iterations: %s' % m.niter)
        print('Fitted pars:%s' % m.params)
        print('Uncertainties: %s' % m.perror)  # TODO: We can use them we define a realistic error on the flux
        print('Value of the summed squared residuals: %s' % m.fnorm)
        print('Number of calls to the function: %s' % m.nfev)
        print('Calculations finished in %s sec' % int(end_time))
        x_s, y_s = func(m.params, atmtype=model, driver='synth', r=r, fout=fout, **kwargs)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_final = sl(x_obs)
        save_synth_spec(x_obs, flux_final, fname='final.spec')
        chi = ((x_obs - flux_final)**2)
        chi2 = np.sum(chi)
        print('This is your chi2 value: '), chi2

        return m.params, x_s, flux_final


def minimize_synth(p0, x_obs, y_obs, r, fout, **kwargs):
    '''Minimize a synthetic spectrum to an observed

     Input
     -----
     p0 : list
       Initial parameters (teff, logg, feh, vt)
     x_obs : ndarray
       Observed wavelength
     y_obs : ndarray
       Observed flux
     r : ndarray
       ranges of the intervals
     fout : str
       Input line list files

     Output
     -----
     params : list
       Final parameters
     x_final : ndarray
       Final wavelength
     y_final : ndarray
       Final synthetic flux
    '''

    from utils import fun_moog_synth as func
    from mpfit import mpfit
    from scipy.interpolate import InterpolatedUnivariateSpline
    from synthetic import save_synth_spec

    model = kwargs['model']
    fix_teff = 1 if kwargs['fix_teff'] else 0
    fix_logg = 1 if kwargs['fix_logg'] else 0
    fix_feh = 1 if kwargs['fix_feh'] else 0
    fix_vt = 1 if kwargs['fix_vt'] else 0
    fix_vmac = 1 if kwargs['fix_vmac'] else 0
    fix_vsini = 1 if kwargs['fix_vsini'] else 0

    def bounds(i, p, model):
        '''Smart way to calculate the bounds of each of parameters'''
        if model.lower() == 'kurucz95':
            bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99, 0, 50, 0, 100]
        if model.lower() == 'apogee_kurucz':
            bounds = [3500, 30000, 0.0, 5.0, -5, 1.5, 0, 9.99, 0, 50, 0, 100]
        if model.lower() == 'marcs':
            bounds = [2500, 8000, 0.0, 5.0, -5, 1.0, 0, 9.99, 0, 50, 0, 100]

        if p[int((i-1)/2)] < bounds[i-1]:
            p[int((i-1)/2)] = bounds[i-1]
        elif p[int((i-1)/2)] > bounds[i]:
            p[int((i-1)/2)] = bounds[i]
        return p

    # Set PARINFO structure for all 6 free parameters for mpfit
    # Teff, logg, feh, vt, vmac, vsini
    # The limits are also cheched by the bounds function
    teff_info  = {'limited': [1, 1], 'limits': [2800.0, 7500.0], 'step': 30,   'mpside': 2, 'fixed': fix_teff}
    logg_info  = {'limited': [1, 1], 'limits': [0.5, 5.0],       'step': 0.2,  'mpside': 2, 'fixed': fix_logg}
    feh_info   = {'limited': [1, 1], 'limits': [-5.0, 1.0],      'step': 0.05, 'mpside': 2, 'fixed': fix_feh}
    vt_info    = {'limited': [1, 1], 'limits': [0.0, 9.99],      'step': 0.3,  'mpside': 2, 'fixed': fix_vt}
    vmac_info  = {'limited': [1, 1], 'limits': [0.0, 50.0],      'step': 0.5,  'mpside': 2, 'fixed': fix_vmac}
    vsini_info = {'limited': [1, 1], 'limits': [0.0, 100.0],     'step': 0.5,  'mpside': 2, 'fixed': fix_vsini}

    parinfo = [teff_info, logg_info, feh_info, vt_info, vmac_info, vsini_info]

    def myfunct(p, x_obs=None, r=None, fout=None, model=None,
                y=None, **kwargs):
        '''Function that return the weighted deviates (to be minimized).

        Input
        ----
        p : list
          Parameters for the model atmosphere
        x_obs : ndarray
          Wavelength
        r : ndarray
          ranges of the intervals
        fout : str
          Line list file
        model : str
          Model atmosphere type
        y : ndarray
          Observed flux


        Output
        -----
        (y-ymodel)/err : ndarray
          Model deviation from observation
        '''

        # Definition of the Model spectrum to be iterated
        options = kwargs['options']
        # Check for bounds
        p = bounds(1, p, model)
        p = bounds(3, p, model)
        p = bounds(5, p, model)
        p = bounds(7, p, model)
        p = bounds(9, p, model)
        p = bounds(11, p, model)
        x_s, y_s = func(p, atmtype=model, driver='synth', r=r, fout=fout, **options)
        sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
        flux_s = sl(x_obs)
        ymodel = flux_s
        # Error on the flux #needs corrections
        err = np.zeros(len(y)) + 0.01
        status = 0
        return([status, (y-ymodel)/err])

    # A dictionary which contains the parameters to be passed to the
    # user-supplied function specified by myfunct via the standard Python
    # keyword dictionary mechanism. This is the way you can pass additional
    # data to your user-supplied function without using global variables.
    fa = {'x_obs': x_obs, 'r': r, 'fout': fout, 'model': model, 'y': y_obs,
          'options': kwargs}

    # Minimization starts here
    # Measure time
    start_time = time.time()
    m = mpfit(myfunct, xall=p0, parinfo=parinfo, ftol=1e-5, xtol=1e-5, gtol=1e-10, functkw=fa)
    print('status = %s' % m.status)
    print('Iterations: %s' % m.niter)
    print('Fitted pars:%s' % m.params)
    print('Uncertainties: %s' % m.perror)  # TODO: We can use them we define a realistic error on the flux
    print('Value of the summed squared residuals: %s' % m.fnorm)
    print('Number of calls to the function: %s' % m.nfev)
    end_time = time.time()-start_time
    print('Calculations finished in %s sec' % int(end_time))
    x_s, y_s = func(m.params, atmtype=model, driver='synth', r=r, fout=fout, **kwargs)
    sl = InterpolatedUnivariateSpline(x_s, y_s, k=1)
    flux_final = sl(x_obs)
    chi = ((x_obs - flux_final)**2)
    chi2 = np.sum(chi)
    print('This is your chi2 value: '), chi2
    # TODO create a header with the parameters in the output file
    save_synth_spec(x_obs, flux_final, fname='final.spec')
    return m.params, x_obs, flux_final


def mcmc_synth(x0, observed, limits):
    '''This could be cool if it worked'''

    import emcee
    from utils import interpol_synthetic, fun_moog as func

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
