#!/usr/bin/env python
# -*- coding: utf8 -*-

# My imports
from __future__ import division, print_function
import logging
import os
from shutil import copyfile
from glob import glob
import numpy as np
import collections

def _run_ares():
    """Run ARES"""
    os.system('ARES > /dev/null')
    for tmp in ['tmp', 'tmp2', 'tmp3']:
        if os.path.isfile(tmp):
            os.remove(tmp)


def round_up0(i):
    """Round up to 2nd decimal because of stupid python"""
    import decimal

    rounded = decimal.Decimal(str(i)).quantize(decimal.Decimal('1.11'), rounding=decimal.ROUND_HALF_UP)
    return np.float(rounded)


def make_linelist(line_file, ares, cut):
    """Merging linelist with ares file
    """
    import pandas as pd

    linelist = pd.read_csv(line_file,  skiprows=2, names=['WL', 'num', 'EP', 'loggf', 'element', 'EWsun'], delimiter=r'\s+')
    linelist = linelist.sort_values(['WL'])

    data = pd.read_csv(ares, usecols=(0,4), names=['wave', 'EW'],  delimiter=r'\s+')
    data = data.sort_values(['wave'])

    dout = pd.merge(left=linelist, right=data, left_on='WL', right_on='wave', how='inner')
    print(dout)
    dout = dout.loc[:, ['WL', 'num', 'EP', 'loggf', 'EW']]
    dout = dout[dout.EW < float(cut)]
    dout.drop_duplicates('WL', keep='first', inplace=True)
    dout = dout.sort_values(['num', 'WL'])

    N2 = len(dout)
    N = len(linelist)
    print('\tARES measured %i lines.' % len(data))
    if N-N2:
        print('\tARES did not find %i lines out of %i of the input line list' % (N-N2, N))
    print('\tThe final line list has %i lines' % N2)

    fout = '%s.moog' % ares.rpartition('.')[0]
    np.savetxt(fout, dout.values, fmt=('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f'), header=' %s' % fout)
    os.remove(ares)
    return


def _options(options=None):
    '''Reads the options inside the config file'''
    defaults = {'lambdai':'3900.0',
                'lambdaf': '25000.0',
                'smoothder': '4',
                'space': '2.0',
                'rejt': False,
                'lineresol': '0.07',
                'miniline': '2.0',
                'plots_flag': False,
                'EWcut': '200.0',
                'snr': False,
                'output': False,
                'rvmask': '"0,0"',
                'force': False
                }
    if not options:
        defaults['rejt'] = '3;5764,5766,6047,6053,6068,6076'
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
        defaults['lineresol'] = float(defaults['lineresol'])
        defaults['miniline'] = float(defaults['miniline'])
        defaults['rvmask'] = str(defaults['rvmask'])
        defaults['EWcut'] = float(defaults['EWcut'])
        if defaults['plots_flag']:
            defaults['plots_flag'] = '1'

        if not defaults['rejt']:
            defaults['rejt'] = '3;5764,5766,6047,6053,6068,6076'
        try:
            defaults['rejt'] = float(defaults['rejt'])
        except ValueError:
            defaults['rejt'] = str(defaults['rejt'])
        return defaults


def update_ares(line_list, spectrum, out, options, fullpath):
    """Driver for ARES"""

    default_options = options
    for key, value in default_options.iteritems():
        if key not in options.keys():
            options[key] = value

    def rejt_from_snr(snr):
        """Calculate rejt from SNR"""
        return 1.0-(1.0/float(snr))

    if options['snr']:
        rejt = rejt_from_snr(options['snr'])
    else:
        rejt = options['rejt']
    plot = 1 if options['plots_flag'] else 0
    if isinstance(rejt, float):
        rejt = 0.999 if rejt > 0.999 else rejt

    if options['output']:
        out = options['output']
    else:
        out = '%s.ares' % spectrum.rpartition('.')[0]

    if fullpath:
        fout = 'specfits=\'%s\'\n' % spectrum
    else:
        fout = 'specfits=\'spectra/%s\'\n' % spectrum
    fout += 'readlinedat=\'rawLinelist/%s\'\n' % line_list
    fout += 'fileout=\'linelist/%s\'\n' % out
    fout += 'lambdai=%s\n' % options['lambdai']
    fout += 'lambdaf=%s\n' % options['lambdaf']
    fout += 'smoothder=%s\n' % options['smoothder']
    fout += 'space=%s\n' % options['space']
    fout += 'rejt=%s\n' % rejt
    fout += 'lineresol=%s\n' % options['lineresol']
    fout += 'miniline=%s\n' % options['miniline']
    fout += 'plots_flag=%s\n' % plot
    fout += 'rvmask=\'0,%s\'\n' % options['rvmask']

    with open('mine.opt', 'w') as f:
        f.writelines(fout)


def findBadLine():
    """Read logARES.txt and return the last measured line (the bad one)"""
    linePresent = False
    with open('logARES.txt', 'r') as lines:
        for line in lines:
            if line.startswith('line result'):
                line = line.split(':')
                badLine = float(line[-1])
                linePresent = True
    if linePresent:
        return badLine
    else:
        return None


def cleanLineList(linelist, badline):
    badline = str(round(badline, 2))
    with open(linelist, 'r') as lines:
        fout = ''
        for line in lines:
            if line.startswith(badline):
                continue
            fout += line
    with open(linelist, 'w') as f:
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

    # Check if there is a directory called linelist, if not create it and ask the user to put files there
    if not os.path.isdir('linelist'):
        os.mkdir('linelist')
        logger.info('linelist directory was created')

    if not os.path.isdir('rawLinelist'):
        os.mkdir('rawLinelist')
        logger.info('linelist directory was created')
        raise IOError('Please put linelists in rawLinelist folder')

    with open(starLines, 'r') as lines:
        for line in lines:
            if not line[0].isalpha():
                logger.debug('Skipping header: %s' % line.strip())
                continue
            logger.info('Processing: %s' % line.strip())
            line = line.strip()
            line = line.split(' ')

            if len(line) == 2:
                options = _options()
                line_list = line[0]
                spectrum = line[1]
            elif len(line) == 3:
                options = _options(line[-1])
                line_list = line[0]
                spectrum = line[1]
            else:
                logger.error('Could not process information for this line: %s' % line)
                continue

            if options['output']:
                out = options['output']
            else:
                out = '%s.ares' % spectrum.rpartition('/')[2].rpartition('.')[0]
                options['output'] = out
            if os.path.isfile('spectra/%s' % spectrum):
                update_ares(line_list, spectrum, out, options, fullpath=False)
            elif os.path.isfile(spectrum):
                 update_ares(line_list, spectrum, out, options, fullpath=True)
            else:
                 continue

            print('Using linelist: %s' % line_list)
            print('Using spectrum: %s' % spectrum)
            if options['output']:
                out = options['output']
                print('Your ARES output: linelist/%s' % out.replace('.ares', '.moog'))
            else:
                out = '%s.ares' % spectrum.rpartition('.')[0]
                print('Your ARES output: linelist/%s' % out.replace('.ares', '.moog'))

            if options['force']:
                index = 1
                while True:
                    _run_ares()
                    if os.path.isfile('linelist/'+out) or os.path.isfile('linelist/'+options['output']):
                        break
                    else:
                        atomicLine = findBadLine()
                        if atomicLine:
                            print('\tRemoving line: %.2f' % atomicLine)
                            copyfile('rawLinelist/'+line_list, 'rawLinelist/tmp%i' % index)
                            line_list = 'tmp%i' % index
                            cleanLineList('rawLinelist/'+line_list, atomicLine)
                            update_ares(line_list, spectrum, out, options)
                            index += 1
                        else:
                            break
                for tmp in glob('rawLinelist/tmp*'):
                    os.remove(tmp)
            else:
                _run_ares()
            line_list = line[0]
            try:
                if options['output']:
                    out = options['output']
                else:
                    out = '%s.ares' % spectrum.rpartition('.')[0]
                make_linelist('rawLinelist/'+line_list, 'linelist/'+out, cut=options['EWcut'])
            except IOError:
                raise IOError('ARES did not run properly. Take a look at "logARES.txt" for more help.')
            print('\n')

    os.remove('logARES.txt')


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        cfgfile = sys.argv[1]
    else:
        cfgfile = 'StarMe_ares.cfg'
    aresdriver(starLines=cfgfile)
