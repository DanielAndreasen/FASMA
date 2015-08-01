#!/usr/bin/python
import argparse

parser = argparse.ArgumentParser(description='This function calculates mic and mac velocities.')
parser.add_argument('teff', help='Input teff')
parser.add_argument('logg', help='Input logg to define evolutionary stage.')
args = parser.parse_args()

if float(args.logg) >= 3.95:   #Dwarfs
    teff = float(args.teff)
    logg = float(args.logg)

    mic = 6.932*teff*(10**-4)-0.348*logg-1.437
    mac = 3.21 + (2.33*(10**-3)*(teff-5777)) + (2.0*(10**-6)*(teff-5777)*(teff-5777)) - (2.0*(logg-4.44)) #Doyle 2014
    print 'mic ', mic
    print 'mac ', mac

if float(args.logg) < 3.95:  #Giants
    teff = float(args.teff)
    logg = float(args.logg)

    mic=3.7-(5.1*teff*(10**-4))
    mac = 3.21 + (2.33*(10**-3)*(teff-5777)) + (2.0*(10**-6)*(teff-5777)*(teff-5777)) - (2.0*(logg-4.44)) #Doyle 2014
    print 'mic ', mic
    print 'mac ', mac


