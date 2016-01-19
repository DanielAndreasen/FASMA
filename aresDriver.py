#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
import yaml
import numpy as np

def _run_ares():
    """Run ARES"""
    os.system('ARES2')


def make_linelist(line_file, ares, cut=200):
    """This function creates a MOOG readable file from the ARES output using the
    line list and atomic data of make_linelist.dat file"""

    # Read the line list and check for multiple identical lines
    linelist = np.genfromtxt(line_file, dtype=None, skiprows=2, names=['line', 'atomic', 'excitation', 'loggf', 'ew_sun'])
    linelist_wave = linelist['line']
    linelist_excitation = linelist['excitation']
    linelist_loggf = linelist['loggf']
    linelist_atomic = linelist['atomic']
    assert (len(np.unique(linelist_wave)) == len(linelist_wave)), 'Check for multiple lines in make_linelist.dat'

    # Read the lines and ews in the ares data and check for identical lines
    data = np.genfromtxt(ares, skip_footer=1, dtype=None, names=['wave', 'fit', 'c1', 'c2', 'ew', 'dew', 'c3', 'c4', 'wave_cal'])
    wave_ares = data['wave']
    ew_ares = data['ew']
    dew_ares = data['dew']
    error_ew = (dew_ares*100)/ew_ares

    assert (len(np.unique(wave_ares)) == len(wave_ares)), 'Check for multiple lines in line.star.ares'

    # Wavelength and EW taken from the ares file.
    # Test whether each element of a 1D array is also present in a second array
    index_ares = np.in1d(wave_ares, linelist_wave)
    common_wave = wave_ares[index_ares]
    ew = ew_ares[index_ares]
    dew = dew_ares[index_ares]
    error_ew = error_ew[index_ares]

    # Sort common elements from ares by wavelength
    ares_values = np.column_stack((common_wave, ew, dew, error_ew))
    ares_sorted = sorted(ares_values, key=lambda row: row[0])
    ares_sorted = np.transpose(ares_sorted)
    indices_1 = ares_sorted[:, 1] > cut
    print('%s lines with EW higher than %s were deleted' % (len(ares_sorted[indices_1]), cut))
    ares_sorted = ares_sorted[~indices_1]
    wave_ares = ares_sorted[0]
    ew = ares_sorted[1]

    # Wavelength and atomic data taken from the make_linelist.dat file.
    # Test whether each element of a 1D array is also present in a second array
    linelist_index = np.in1d(linelist_wave, wave_ares)
    index_lines_not_found = np.invert(np.in1d(linelist_wave, wave_ares))
    lines_not_found = linelist_wave[index_lines_not_found]
    print('ARES did not find ', len(lines_not_found), 'lines: ', lines_not_found)
    wave = linelist_wave[linelist_index]
    excitation = linelist_excitation[linelist_index]
    loggf = linelist_loggf[linelist_index]
    atomic = linelist_atomic[linelist_index]
    print('Lines in the new line list: ', len(wave))
    # Sort common elements from line list by wavelength
    linelist_values = np.column_stack((wave, atomic, excitation, loggf))
    linelist_sorted = sorted(linelist_values, key=lambda row: row[0])
    linelist_sorted = np.transpose(linelist_sorted)

    # Merge line list data with the EW from ARES
    # Sort the FeI and the FeII lines using the atomic number
    values = np.column_stack((ares_sorted[0], linelist_sorted[1], linelist_sorted[2], linelist_sorted[3], ares_sorted[1]))
    sorted_values = sorted(values, key=lambda row: row[1])
    sorted_values = np.transpose(sorted_values)

    # Write results in MOOG readable format
    assert np.array_equal(ares_sorted[0], linelist_sorted[0]), 'There is something wrong with the common elements of ARES and the line list'
    data = zip(sorted_values[0], sorted_values[1], sorted_values[2], sorted_values[3], sorted_values[4])
    np.savetxt('%s.moog' % ares, data, fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f'), header=' %s' % ares)


def _options(options=False):
    '''Reads the options inside the config file'''
    defaults = {'lambdai':'3900.0',
                'fileout': False,
                'lambdaf': '25000.0',
                'smoothder': '4',
                'space': '2.0',
                'rejt': '0.995',
                'lineresol' : '0.07',
                'miniline' : '2.0',
                'plots_flag': '0',
                'EWcut': '200.0',
                'snr': False,
                'rvmask' : '"0,0"'
                }
    if not options:
        return defaults
    else:
        options = options.split(',')
        for option in options:
            if ':' in option:
                option = option.split(':')
                defaults[option[0]] = option[1]
            else:
                defaults[option] = True
        defaults['lambdai'] = float(defaults['lambdai'])
        defaults['lambdaf'] = float(defaults['lambdaf'])
        defaults['smoothder'] = int(defaults['smoothder'])
        defaults['space'] = float(defaults['space'])
        defaults['rejt'] = float(defaults['rejt'])
        defaults['lineresol'] = float(defaults['lineresol'])
        defaults['miniline'] = float(defaults['miniline'])
        defaults['plots_flag'] = int(defaults['plots_flag'])
        defaults['rvmask'] = str(defaults['rvmask'])
        defaults['EWcut'] = float(defaults['EWcut'])
        return defaults


def update_ares(line_list, spectrum, out, options):
    """Driver for ARES"""

    default_options = options
    for key, value in default_options.iteritems():
        if key not in options.keys():
            options[key] = value

    def rejt_from_snr(snr):
        """Calculate rejt from SNR"""
        return 1.0-(1.0/snr)

    if options['snr']:
        rejt = rejt_from_snr(options['snr'])
    else:
        rejt = options['rejt']
    rejt = 0.999 if rejt > 0.999 else rejt
    plot = 1 if options['plots_flag'] else 0

    fout = 'specfits=\'spectra/%s\'\n' % spectrum
    fout += 'readlinedat=\'linelist/%s\'\n' % line_list
    fout += 'fileout=\'linelist/%s\'\n' % out
    fout += 'lambdai=%s\n' % options['lambdai']
    fout += 'lambdaf=%s\n' % options['lambdaf']
    fout += 'smoothder=%s\n' % options['smoothder']
    fout += 'space=%s\n' % options['space']
    fout += 'rejt=%s\n' % rejt
    fout += 'lineresol=%s\n' % options['lineresol']
    fout += 'miniline=%s\n' % options['miniline']
    fout += 'plots_flag=%s\n' % plot

    with open('mine.opt', 'w') as f:
        f.writelines(fout)


def aresdriver(starLines='StarMe_ares.cfg'):
    """The function that glues everything together

    Input:
    starLines   -   Configuration file (default: StarMe_ares.cfg)

    Output:
    <linelist>.out          -   Output file
    """
    try:  # Cleaning from previous runs
        os.remove('captain.log')
    except OSError:
        pass
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler('captain.log')
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    #Check if there is a directory called linelist, if not create it and ask the user to put files there
    if not os.path.isdir('linelist'):
        os.mkdir('linelist')
        logger.info('linelist directory was created')

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Line list: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            #Check if the linelist is inside the directory if not log it and pass to next linelist
            if not os.path.isfile(line[0]):
                logger.error('Error: %s not found.' % line[0])
                continue

            if len(line) == 3:
                options = _options()
                line_list = line[0]
                spectrum = line[1]
                out = line[2]
                update_ares(line_list, spectrum, out, options)
            elif len(line) == 4:
                line_list, spectrum, out = map(str, line[0:-1])
                options = _options(line[-1])
                update_ares(line_list, spectrum, out, options)
            else:
                logger.error('Could not process information for this line: %s' % line)
                continue

            _run_ares()
            line_list = line[0]
            out = 'linelist/%s' % line[2]
            make_linelist(line_list, out, cut=options['EWcut'])

if __name__ == '__main__':
    aresdriver(starLines='StarMe_ares.cfg')
