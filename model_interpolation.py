#!/usr/bin/python
import numpy as np
from numpy import *
from scipy.integrate import simps
import scipy as sp


"""
Following the concept of Sz. Mezeros and C. Allende Prieto:
For each set of parameters, we identified the 8 immediate neighbors with higher and lower values for each parameter in
the grid, calculated by numerical integration the Rosseland optical depth for each, re-sampled all the
thermodynamical quantities in the atmosphere (temperature, gas pressure, and electron density) on
a common optical depth scale for all models by linear interpolation, and then interpolated, linearly,
all the thermodynamical quantities to the parameters (Teff , log g, and [Fe/H]) of the target model.
Other quantities included in the models (Rosseland opacities, radiative pressure, etc.) were
also interpolated in the same way.
"""
#def t_ross_scale(model_rhox, model_abross, num_layers):
#       tauross = [] 
#	for i in range(int(num_layers)):
#	        tauross.append(sp.integrate.simps(model_rhox[i], model_abross[i]))
#        tauross = np.array(tauross)
#        return tauross

def read_model(filename):

	with open(filename) as f:
    #First I read the header 
		lines = f.readlines()[:23]
	for line in lines:		
		vline = line.split()
		param = vline[0]
		if param == "TEFF":
			teff = float(vline[1])
			logg = float(vline[3])
			#print teff, logg
    elif param == "READ":
      num_layers = int(vline[2])
			#print num_layers
  #This are the thermodynamical quantities.
	model = np.genfromtxt(filename, dtype=None, skiprows=23, skip_footer=2, names = ['RHOX', 'T', 'P', 'XNE', 'ABROSS', 'ACCRAD', 'VTURB', 'tab1', 'tab2'])
	model_rhox = model['RHOX']
	model_t = model['T']
	model_p = model['P']
	model_xne = model['XNE']
	model_abross = model['ABROSS']
	model_accrad = model['ACCRAD']
	model_vturb = model['VTURB']
	model_tab1 = model['tab1']
	model_tab2 = model['tab2']

	if len(model) != num_layers:
	        raise Exception("FORMAT ERROR")
  #I want to build tauross.. But I dont know exactly what I am doing..
  tauross = []
  tauross_0 = model_abross[0]*model_rhox[0] #This is supposed to be the first element
  for i in range(1,num_layers):
       	tauross.append(sp.integrate.simps(model_rhox[0:i], model_abross[0:i],  even='last'))
  tauross = np.array(tauross)
  tauross = np.insert(tauross, 0, tauross_0)
	#tauross = t_ross_scale(model_rhox, model_abross, num_layers)	

  return model_rhox, model_t, model_p, model_xne, model_abross, model_accrad, model_vturb, tauross

#x1 = read_model('5500g45.p01')
#x2 = read_model('5750g45.p01')
#x3 = read_model('6000g45.p01')
#x4 = read_model('6250g45.p01')
#print x1[2]
#print x2[2]
#print x3[2]
#print x4[2]

avail_models = ['5500g45.p01', '5750g45.p01', '6000g45.p01', '6250g45.p01'] 
num_of_models = len(avail_models)

#Now I will resample everything on the rosseland optical depth scale of model1 

tauross_max = []
tauross_min = []

for i, x in enumerate(avail_models):
     models.append(read_model(x))
     tauross_min.append(read_model(x)[7][0])
     tauross_max.append(read_model(x)[7][-1])

max_tauross = max(tauross_max)
min_tauross = min(tauross_min)
print max_tauross
