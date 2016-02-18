#!/usr/bin/python

import numpy as np
import pandas as pd
import seaborn as sns
sns.set_style('white')
sns.set_context('talk', font_scale=1.2)
import matplotlib.pyplot as plt

def massTorres(teff, erteff, logg, erlogg, feh, erfeh):
    """Calculate a mass using the Torres calibration"""
    ntrials = 10000
    randomteff = teff + erteff * np.random.randn(ntrials)
    randomlogg = logg + erlogg * np.random.randn(ntrials)
    randomfeh = feh + erfeh * np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    a1 = 1.5689
    a2 = 1.3787
    a3 = 0.4243
    a4 = 1.139
    a5 = -0.1425
    a6 = 0.01969
    a7 = 0.1010

    M = np.zeros(ntrials)
    logM = np.zeros(ntrials)
    for i in xrange(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logMass = a1 + a2 * X + a3 * X * X + a4 * X * X * X + a5 *\
            randomlogg[i] * randomlogg[i] + a6 * randomlogg[i] *\
            randomlogg[i] * randomlogg[i] + a7 * randomfeh[i]
        logM[i] = logMass
        M[i] = 10 ** logMass

    meanMasslog = np.mean(logM)
    sigMasslog = np.sqrt(np.sum([(logMi - meanMasslog)**2 for logMi in logM]) /
                         (ntrials - 1))
    sigMasslogTot = np.sqrt(0.027*0.027 + sigMasslog*sigMasslog)

    meanMass = 10**meanMasslog
    sigMass = 10**(meanMasslog + sigMasslogTot) - meanMass

    return meanMass, sigMass


def radTorres(teff, erteff, logg, erlogg, feh, erfeh):
    ntrials = 10000
    randomteff = teff + erteff*np.random.randn(ntrials)
    randomlogg = logg + erlogg*np.random.randn(ntrials)
    randomfeh = feh + erfeh*np.random.randn(ntrials)

    # Parameters for the Torres calibration:
    b1 = 2.4427
    b2 = 0.6679
    b3 = 0.1771
    b4 = 0.705
    b5 = -0.21415
    b6 = 0.02306
    b7 = 0.04173

    R = np.zeros(ntrials)
    logR = np.zeros(ntrials)

    for i in xrange(ntrials):
        X = np.log10(randomteff[i]) - 4.1
        logRad = b1 + b2 * X + b3 * X * X + b4 * X * X * X + b5 *\
            randomlogg[i] * randomlogg[i] + b6 * randomlogg[i] *\
            randomlogg[i] * randomlogg[i] + b7 * randomfeh[i]
        logR[i] = logRad
        R[i] = 10 ** logRad

    meanRadlog = np.mean(logR)
    sigRadlog = np.sqrt(np.sum([(logRi-meanRadlog)**2 for logRi in logR]) /
                        (ntrials-1))
    sigRadlogTot = np.sqrt(0.014*0.014 + sigRadlog*sigRadlog)

    meanRad = 10**meanRadlog
    sigRad = 10**(meanRadlog + sigRadlogTot) - meanRad

    return meanRad, sigRad


if __name__ == '__main__':
    df = pd.read_csv('results.csv', delimiter=r'\s+')
    m = [massTorres(t, et, l, el, f, ef) for t, et, l, el, f, ef in zip(df.teff, df.tefferr, df.logg, df.loggerr, df.feh, df.feherr)]
    r = [radTorres(t, et, l, el, f, ef) for t, et, l, el, f, ef in zip(df.teff, df.tefferr, df.logg, df.loggerr, df.feh, df.feherr)]
    df['mass'] = pd.Series(np.asarray(m)[:, 0])
    df['masserr'] = pd.Series(np.asarray(m)[:, 1])
    df['radius'] = pd.Series(np.asarray(r)[:, 0])
    df['radiuserr'] = pd.Series(np.asarray(r)[:, 1])
    df['lum'] = (df.teff/5777)**4 * df.mass

    # Plot the results
    color = sns.color_palette()
    plt.figure()
    plt.plot(df.teff[df.convergence], df.lum[df.convergence], 'o', c=color[0], ms=10)
    plt.plot(df.teff[~df.convergence], df.lum[~df.convergence], 'o', c=color[2], ms=5)
    plt.xlim(plt.xlim()[::-1])  # Reverse Teff (as normal HR diagram)
    plt.title('HR diagram')
    plt.xlabel(r'Teff [K]')
    plt.ylabel(r'L$_\odot$')
    plt.tight_layout()
    plt.show()
