#!/usr/bin/python

import numpy as np
import urllib


def mass_radius_padova(star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal):
	"""Returns mass and radius from padova interface.
	   Enter star, vmag, parallax, er_parallax, temp, er_temp, metal, er_metal
	"""
	url = 'http://stev.oapd.inaf.it/cgi-bin/param'
	#These are the parameters in the webpage to tune
	form_data = {'param_version': '1.3',
	             'star_name': star,
				 'star_teff': temp,
				 'star_sigteff': er_temp,
				 'star_feh': metal,
				 'star_sigfeh': er_metal,
				 'star_vmag': vmag,
				 'star_sigvmag': '0.02',
				 'star_parallax': parallax,
				 'star_sigparallax': er_parallax,
				 'isoc_kind': 'parsec_CAF09_v1.1',
				 'kind_interp': '1',
				 'kind_tpagb': '0',
				 'kind_pulsecycle': '0',
				 'kind_postagb': '-1',
				 'imf_file': 'tab_imf/imf_chabrier_lognormal.dat',
				 'sfr_file': 'tab_sfr/sfr_const_z008.dat',
				 'sfr_minage': '0.1e9',
				 'sfr_maxage': '12.0e9',
				 'flag_sismo': '0',
				 'submit_form': 'Submit'}
	urllib.urlretrieve(url, 'parameters.html', lambda x,y,z:0, urllib.urlencode(form_data))

	#write results
	with open('parameters.html') as f:
        line = f.readlines()[19]

	mass = []
	radius = []
	name = []
	ermass = []
	erradius = []

    line = line.replace(' .<p>Results for ', '')
 	line = line.replace('Mass=', '')
 	line = line.replace('&plusmn;', ' , ')
	line = line.replace(' ','')
	line = line.replace('<i>R</i>=','')
	line = line.replace(':Age=',',')
	line = line.replace('<i>M</i>&#9737','')
	line = line.replace('<i>R</i>&#9737','')
 	line = line.split(',')
    name.append(str(line[0]))
    mass.append(str(line[3]))
   	ermass.append(str(line[4]))
   	radius.append(str(line[7]))
   	erradius.append(str(line[8]))

	return name, mass, ermass, radius, erradius


def mass_radius_torres(T, logg, fe, dT, dlogg, dfe):
       """Derive mass and radius using the calibration of Torres 2010 and
	      correction for mass deviations according to Santos et al. 2013"""
    a1 = 1.5689
	a2 = 1.3787
	a3 = 0.4243
	a4 = 1.139
	a5 = -0.1425
	a6 = 0.01969
	a7 = 0.1010
	a = [1.5689, 1.3787, 0.4243, 1.139, -0.1425, 0.01969, 0.1010]

	b1 = 2.4427
	b2 = 0.6679
	b3 =  0.1771
	b4 = 0.705
	b5 = -0.21415
	b6 = 0.02306
	b7 = 0.04173
	da1 = 0.058
	da2 = 0.029
	da3 = 0.029
	da4 = 0.240
	da5 = 0.011
	da6 = 0.0019
	da7 = 0.014
	db1 = 0.038
	db2 = 0.016
	db3 = 0.027
	db4 = 0.130
	db5 = 0.0075
	db6 = 0.0013
	db7 = 0.0082

    X = np.log10(T)-4.1
    dX = dT/T
	log_M = np.polyval(a[0:4], X) + a[4]*logg**2 + a[5]*logg**3 + a[6]*fe
	dlog_M = np.sqrt((da1**2) + ((a2*dX)**2) + ((da2*X)**2) + ((da3*(X**2))**2) + ((a3*2*X*dX)**2) + ((da4*(X**3))**2) + ((a4*3*X*X*dX)**2) + ((da5*(logg**2))**2) + ((a5*2*logg*dlogg)**2) + ((da6*(logg**3))**2) + ((a6*3*logg*logg*dlogg)**2) + ((da7*fe)**2) + ((a7*dfe)**2))
    log_R = b1 + b2*X + b3*(X**2) + b4*(X**3) + b5*(logg**2) + b6*(logg**3) + b7*fe
    dlog_R = np.sqrt((db1**2) + ((b2*dX)**2) + ((db2*X)**2) + ((db3*(X**2))**2) + ((b3*2*X*dX)**2) + ((db4*(X**3))**2) + ((b4*3*X*X*dX)**2) + ((db5*(logg**2))**2) + ((b5*2*logg*dlogg)**2) + ((db6*(logg**3))**2) + ((b6*3*logg*logg*dlogg)**2) + ((db7*fe)**2) + ((b7*dfe)**2))
    Mt = 10**log_M
    Rt = 10**log_R
    dMt = dlog_M*Mt * 10**(log_M-1.0)
    dRt = dlog_R*Rt * 10**log_R-1.0)
    Mcal = 0.791*Mt*Mt - 0.575*Mt +0.701
    dMcal = np.sqrt(((0.791*Mt*dMt)**2) + ((0.575*dMt)**2))
    return Mt, dMt, Rt, dRt, Mcal, dMcal

def bc_correction(teff):
    """Bolometric correction according to Flower 1996 and Torres 2010"""
    lteff= np.log10(teff)
	c1 = [-0.190537291496456e+05, 0.155144866764412e+05, -0.421278819301717e+04, 0.381476328422343e+03]
	c2 = [-0.370510203809015e+05, 0.385672629965804e+05, -0.150651486316025e+05, 0.261724637119416e+04, -0.170623810323864e+03]
	c3 = [-0.370510203809015e+05, 0.385672629965804e+05, -0.150651486316025e+05, 0.261724637119416e+04, -0.170623810323864e+03]
    if lteff<3.7:
		return np.polyval(c1, teff)
    elif 3.7<=lteff and lteff<3.9:
		return np.polyval(c2, teff)
	elif 3.9<=lteff:
		return np.polyval(c3, teff)

def logg_trigomonetric(teff, mass, v, bc, par, dpar, dteff, dmass):
    """Calculate logg from parallax"""
    logg = 4.44 + np.log10(mass) + 4.*np.log10(teff/5777.) + 0.4*(v+bc) + 2.*np.log10(par/1000.)+0.108
    logg = np.round(logg,2)
    dlogg = np.sqrt((dmass/mass)**2 + (dteff/teff)**2 + 2*(dpar/par)**2)
    dlogg = np.round(dlogg,2)
    return logg, dlogg

def logg_transit_mortier(logg, teff):
    """Correct logg according to Mortier et al. 2014 using transit data
       Valid for 4500 < teff < 6500"""
    return -4.57e-4*teff+2.59+logg

def logg_astero_mortier(logg, teff):
    """Correct logg according to Mortier et al. 2014 using asteroseismic data
    Valid for 5000 < teff < 6500"""
	return -3.89e-4*teff+2.10+logg