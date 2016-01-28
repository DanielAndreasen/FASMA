#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division
import numpy as np
import os
from copy import copy


def print_format(iter, x, slopes):
    """Print the stellar atmospheric parameters in a nice format"""
    print '{:4d}{:>6d}{:>8.2f}{:>+9.2f}{:>8.2f}{:>+9.3f}{:>+11.3f}{:>11.2f}'.format(iter, x[0], x[1], x[2], x[3], slopes[0], slopes[1], slopes[2])


def check_bounds(parameter, bounds, i):
    """
    Function which checks if parameters are within bounds.
    Input - parameter: what we want to check; bounds: ze bounds;
    i: the index of the bounds we want to check"""
    if parameter < bounds[i-1]:
        parameter = bounds[i-1]
    elif parameter > bounds[i]:
        parameter = bounds[i]
    return parameter


def check_convergence(RW, EP, Abdiff, fe_input, fe):
    """
    Check convergence criteria
    """
    EP = 0.00 if f_teff else EP
    RW = 0.00 if f_vt else RW
    Abdiff = 0.00 if f_logg else Abdiff
    fe = fe_input+7.47 if f_feh else fe

    cond1 = abs(RW) <= RWcriteria
    cond2 = abs(Abdiff) <= ABdiffcriteria
    cond3 = abs(EP) <= EPcriteria
    cond4 = round(fe_input, 2) == round(fe-7.47, 2)
    return cond1 and cond2 and cond3 and cond4


def _bump(x, alpha):
    """Bump to the values in the list, x"""
    for i, X in enumerate(zip(alpha, x)):
        ai, xi = X
        sig = 0.01 if ai*xi == 0 else ai*xi
        if ai:
            x[i] = np.random.normal(xi, abs(sig))
    x[0] = int(x[0])
    x[1] = round(x[1], 2)
    x[2] = round(x[2], 2)
    x[3] = abs(round(x[3], 2))  # We can't have negative microturbulence
    return x


def minimize(x0, func, model="kurucz95", weights='null',
            fix_teff=False, fix_logg=False, fix_feh=False, fix_vt=False,
            iterations=160, EPcrit=0.001, RWcrit=0.003, ABdiffcrit=0.01,
            MOOGv=2013):
    """
    Sane minimization like a normal human being would do it.
    """
    global EPcriteria, RWcriteria, ABdiffcriteria
    global f_teff, f_logg, f_vt, f_feh
    EPcriteria, RWcriteria, ABdiffcriteria = EPcrit, RWcrit, ABdiffcrit
    f_teff, f_logg, f_feh, f_vt = fix_teff, fix_logg, fix_feh, fix_vt
    # Step size in Teff, logg, vt
    step = (700, 1.50, 0.50)
    if model.lower() == "kurucz95":
        bounds = [3750, 39000, 0.0, 5.0, -3, 1, 0, 9.99]

    res, slopeEP, slopeRW, abundances = func(x0, version=MOOGv)
    Abdiff = np.diff(abundances)[0]
    parameters = list(x0)
    if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
        return parameters, True
    try:
        os.remove('minimization_profile.dat')
    except OSError:
        pass

    all_params = [copy(parameters)]
    N = 0
    print(' i    Teff    logg    [Fe/H]    vt    EPslope    RWslope    |FeI-FeII|')
    print('-' * 70)
    while N < iterations:
        # Step for Teff
        if (abs(slopeEP) >= EPcriteria) and not fix_teff:
            s = np.sign(slopeEP)
            step_i = s * step[0]/abs(np.log(abs(slopeEP)+0.0005))**3
            if abs(step_i) < 1:
                step_i = s*1  # Minimum 1K
            elif abs(step_i) > 17625:
                step_i = s*17625  # Maximum 17625K (half the interval?)
            parameters[0] += step_i
            parameters[0] = check_bounds(parameters[0], bounds, 1)
            parameters[0] = int(parameters[0])

        # Step for VT
        if (abs(slopeRW) >= RWcriteria) and not fix_vt:
            s = np.sign(slopeRW)
            step_i = s * step[2]/abs(np.log(abs(slopeRW)+0.0005))**3
            if abs(step_i) < 0.01:
                step_i = s*0.01  # Minimum vt step
            elif abs(step_i) > 2.5:
                step_i = s*2.5  # Maximum vt step
            else:
                step_i = step_i
            parameters[3] += step_i
            parameters[3] = 0.00 if parameters[3] < 0 else parameters[3]
            parameters[3] = check_bounds(parameters[3], bounds, 7)
            parameters[3] = round(parameters[3], 2)

        # Step for logg
        if (abs(Abdiff) >= ABdiffcrit) and not fix_logg:
            s = -np.sign(Abdiff)
            step_i = s * step[1]/abs(np.log(abs(Abdiff)+0.0005))**3
            if abs(step_i) < 0.01:
                step_i = s*0.01  # Minimum logg step
            elif abs(step_i) > 2.5:
                step_i = s*2.5  # Maximum logg step
            parameters[1] += step_i
            parameters[1] = check_bounds(parameters[1], bounds, 3)
            parameters[1] = round(parameters[1], 2)

        # Step for [Fe/H]
        if not fix_feh:
            parameters[2] = abundances[0]-7.47
            parameters[2] = check_bounds(parameters[2], bounds, 5)
            parameters[2] = round(parameters[2], 2)

        if parameters in all_params:
            alpha = [0] * 4
            alpha[0] = abs(slopeEP) if not fix_teff else 0
            alpha[1] = abs(Abdiff) if not fix_logg else 0
            alpha[2] = 0.01 if not fix_feh else 0
            alpha[3] = abs(slopeRW) if not fix_vt else 0
            parameters = _bump(parameters, alpha)
            parameters[0] = int(check_bounds(parameters[0], bounds, 1))
            parameters[1] = check_bounds(parameters[1], bounds, 3)
            parameters[2] = check_bounds(parameters[2], bounds, 5)
            parameters[3] = check_bounds(parameters[3], bounds, 7)
        all_params.append(copy(parameters))

        res, slopeEP, slopeRW, abundances = func(parameters, weights=weights, version=MOOGv)
        Abdiff = np.diff(abundances)[0]
        N += 1
        print_format(N, parameters, (slopeEP, slopeRW, Abdiff))

        if check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0]):
            print 'Stopped in %i iterations' % N
            return parameters, True

    print 'Stopped in %i iterations' % N
    c = check_convergence(slopeRW, slopeEP, Abdiff, parameters[2], abundances[0])
    return parameters, c


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
    import matplotlib.pyplot as plt
    x0 = [5777, 4.44, 0.0, 1.0]
    fname = 'spectra/solar_synthetic.txt'
    limits = (6444.672, 6447.340)

    # MCMC example
    mcmc_synth([5777, 4.44], fname, limits)



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
