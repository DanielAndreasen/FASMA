#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from copy import copy


class Minimize:
    """Minimize for best parameters given a line list"""

    def __init__(self, x0, func, model="kurucz95", weights='null',
                 fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
                 iterations=160, EPcrit=0.001, RWcrit=0.003, ABdiffcrit=0.01,
                 MOOGv=2013, GUI=True):
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


    def _getMic(self):
        """Get the microturbulence if this is fixed"""
        if self.x0[1] >= 3.95:
            self.x0[3] = round(6.932*self.x0[0]/10000 - 0.348*self.x0[1] - 1.437, 2)
        else:
            self.x0[3] = round(3.7 - 5.1*self.x0[0]/10000, 2)


    def print_format(self):
        """Print the stellar atmospheric parameters in a nice format"""
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
        """
        Function which checks if parameters are within bounds.
        Input - parameter: what we want to check; bounds: ze bounds;
        i: the index of the bounds we want to check"""
        if self.x0[int((i-1)/2)] < self.bounds[i-1]:
            self.x0[int((i-1)/2)] = self.bounds[i-1]
        elif self.x0[int((i-1)/2)] > self.bounds[i]:
            self.x0[int((i-1)/2)] = self.bounds[i]
        return self.x0[int((i-1)/2)]


    def check_convergence(self, fe_input):
        """Check the convergence criteria"""
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
        """Bump to the values in the list, x"""
        for i, X in enumerate(zip(alpha, self.x0)):
            ai, xi = X
            sig = 0.01 if ai*xi == 0 else ai*xi
            if ai:
                self.x0[i] = np.random.normal(xi, abs(sig))
        self.x0[0] = int(self.x0[0])
        self.x0[1] = round(self.x0[1], 2)
        self.x0[2] = round(self.x0[2], 2)
        self.x0[3] = abs(round(self.x0[3], 2))


    def minimize(self):
        step = (700, 1.50, 0.50)

        res, self.slopeEP, self.slopeRW, abundances = self.func(self.x0, version=self.MOOGv)
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
                self.x0[0] = self.check_bounds(1)
                self.x0[0] = int(self.x0[0])

            # Step for VT
            if (abs(self.slopeRW) >= self.RWcrit) and not self.fix_vt:
                s = np.sign(self.slopeRW)
                step_i = s * step[2]/abs(np.log(abs(self.slopeRW)+0.0005))**3
                step_i = s*0.01 if abs(step_i) < 0.01 else step_i
                self.x0[3] += step_i
                self.x0[3] = self.check_bounds(7)
                self.x0[3] = round(self.x0[3], 2)

            # Step for logg
            if (abs(self.Abdiff) >= self.ABdiffcrit) and not self.fix_logg:
                s = -np.sign(self.Abdiff)
                step_i = s * step[1]/abs(np.log(abs(self.Abdiff)+0.0005))**3
                step_i = s*0.01 if abs(step_i) < 0.01 else step_i
                self.x0[1] += step_i
                self.x0[1] = self.check_bounds( 3)
                self.x0[1] = round(self.x0[1], 2)

            # Step for [Fe/H]
            if not self.fix_feh:
                self.x0[2] = abundances[0]-7.47
                self.x0[2] = self.check_bounds(5)
                self.x0[2] = round(self.x0[2], 2)

            if self.fix_vt:
                self._getMic()  # Reset the microturbulence
                self.x0[3] = self.check_bounds(7)

            if self.x0 in parameters:
                alpha = [0] * 4
                alpha[0] = abs(self.slopeEP) if not self.fix_teff else 0
                alpha[1] = abs(self.Abdiff) if not self.fix_logg else 0
                alpha[2] = 0.01 if not self.fix_feh else 0
                alpha[3] = abs(self.slopeRW) if not self.fix_vt else 0
                self._bump(alpha)
                self.x0[0] = int(self.check_bounds(1))
                self.x0[1] = self.check_bounds(3)
                self.x0[2] = self.check_bounds(5)
                self.x0[3] = self.check_bounds(7)
            parameters.append(copy(self.x0))

            res, self.slopeEP, self.slopeRW, abundances = self.func(self.x0, weights=self.weights, version=self.MOOGv)
            self.Abdiff = np.diff(abundances)[0]
            self.iteration += 1
            self.print_format()
            best[res] = parameters[-1]
            if self.check_convergence(abundances[0]):
                print 'Stopped in %i iterations' % self.iteration
                return self.x0, True

        print 'Stopped in %i iterations' % self.iteration
        if self.check_convergence(abundances[0]):
            return self.x0, True
        else:
            # Return the best solution rather than the last iteration
            _ = self.func(best[min(best.keys())], weights=self.weights, version=self.MOOGv)
            return best[min(best.keys())], False


def chi2(obs, theory):
    """chi^2 function"""
    error = 1.0
    chi = ((obs - theory)/error)**2
    chi2 = np.sum(chi)
    return chi2


def best_chi(x0, iter_step, steps, wave, flux_obs, limits, p):
    idx = {'teff': 0, 'logg': 1}
    idx = idx[p]
    values = np.linspace(x0[idx]-steps[idx], x0[idx]+steps[idx], 5)

    results = np.zeros((5, 2))
    for value in values:
        if p == 'teff':
            func((value, x0[1], x0[2], x0[3]), driver='synth')
            wavel, flux_obs, flux_inter_synth = interpol_synthetic(wave, flux_obs, limits[0], limits[1])
            chi = chi2(flux_obs, flux_inter_synth)
            results.append((value, chi))
        elif p == 'logg':
            func((x0[0], value, x0[2], x0[3]), driver='synth')
            wavel, flux_obs, flux_inter_synth = interpol_synthetic(wave, flux_obs, limits[0], limits[1])
            chi = chi2(flux_obs, flux_inter_synth)
            results.append((value, chi))
    results = np.array(results)
    chi_best = results[results[:,-1]==results[:,-1].min()]
    return chi_best, results


def minimize_synth(x0, observed, limits):
    '''Minimize a synthetic spectrum to an observed

    Input:
        args:
    Output:
        output
    '''
    from utils import interpol_synthetic, fun_moog as func
    from utils import read_observations

    wavelength_obs, flux_obs = np.loadtxt(observed, unpack=True, usecols=(0, 1))
    idx = (wavelength_obs >= limits[0]) & (wavelength_obs <= limits[1])
    wavelength_obs, flux_obs = wavelength_obs[idx], flux_obs[idx]

    flux_obs /= np.median(flux_obs)
    # Normalization (use first 50 points below 1.2 as constant continuum)
    maxes = flux_obs[(flux_obs < 1.2)].argsort()[-50:][::-1]
    flux_obs /= np.median(flux_obs[maxes])


    # This is dangerous. We have to search all parameter space unless specified by the user.
    func(x0, driver='synth')
    wavelength_obs, flux_obs, flux_inter_synth = interpol_synthetic(wavelength_obs, flux_obs, limits[0], limits[1])
    chi_i = chi2(flux_obs, flux_inter_synth)

    # TODO: ITERATION STARTS, A CONVERGENCE IS NEEDED
    # Maybe the space of search should be defined by the user depending on how well the initial conditions are.
    iterations = 0
    steps = np.array([800.0, 0.9])
    # This is dangerous. We have to search all parameter space unless specified by the user.
    iter_step = steps/(iterations+1)

    tmp = x0[:]
    Teff_rng = np.linspace(x0[0]-50, x0[0]+50, 20)
    logg_rng = np.linspace(x0[1]-0.1, x0[1]+0.1, 10)


    Teffs, loggs = np.meshgrid(Teff_rng, logg_rng)
    z = np.zeros(Teffs.T.shape)
    for i, Teff in enumerate(Teffs[0, :]):
        for j, logg in enumerate(loggs[:, 0]):
            print 'Python --\tTeff: %i/%i\tlogg: %.2f/%.2f' % (Teff, Teff_rng.max(), logg, logg_rng.max())
            func((Teff, logg, x0[2], x0[3]), driver='synth')
            wavel, flux_obs, flux_inter_synth = interpol_synthetic(wavelength_obs, flux_obs, limits[0], limits[1])
            chi = chi2(flux_obs, flux_inter_synth)
            z[i, j] = chi
    return z




    # return d



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







if __name__ == '__main__':
    a = None
    # import matplotlib.pyplot as plt
    # x0 = [5777, 4.44, 0.0, 1.0]
    # fname = 'spectra/solar_synthetic.txt'
    # limits = (6444.672, 6447.340)
    #
    # # MCMC example
    # mcmc_synth([5777, 4.44], fname, limits)



    # z = minimize_synth(x0, fname, limits=limits)
    #
    # z = z.T
    # z1 = z.copy()
    # z /= z.max()
    #
    # Teff_rng = np.linspace(x0[0]-50, x0[0]+50, 20)
    # logg_rng = np.linspace(x0[1]-0.1, x0[1]+0.1, 10)
    #
    #
    # Teffs, loggs = np.meshgrid(Teff_rng, logg_rng)
    # plt.contourf(Teffs, loggs, z, levels=np.linspace(0, 1, 100), cmap=plt.cm.spectral)
    # plt.colorbar()
    # plt.show()
