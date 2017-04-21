import numpy as np
import pandas as pd


def save_loggf(fname, df, linelist):
    '''
    Update the line list with the new loggf values

    Inputs
    ------
    fname : str
      Name of the complete line list in rawLinelist
    df : pandas DataFrame
      The line list loaded as a DataFrame for easy merging
    linelist : str
      The name of the line list, which will be overwritten
    '''

    def merge_loggf(fname, df):
        '''
        Merge the line list DataFrame with the complete linelist

        Inputs
        ------
        fname : str
          Name of the complete line list in rawLinelist
        df : pandas DataFrame
          The line list loaded as a DataFrame for easy merging

        Output
        ------
        df : pandas DataFrame
          The line list with updated loggf values
        '''
        cols = ('wavelength', 'X', 'EP', 'new_gf', 'el', 'ewsun')
        df1 = pd.read_csv(fname, skiprows=2, delimiter=r'\s+',
                          names=cols)
        d = pd.merge(df1, df, left_on='wavelength', right_on='wavelength')
        if len(d) != len(df):
            return None
        else:
            return d.loc[:, ['wavelength', 'X_x', 'EP_x', 'new_gf', 'EW']]

    fname = 'rawLinelist/%s' % fname
    with open(linelist, 'r') as f:
        hdr = f.readline().strip()

    df = merge_loggf(fname, df)
    if df is not None:
        fmt = ('%9.3f', '%10.1f', '%9.2f', '%9.3f', '%28.1f')
        np.savetxt(linelist, df.values, fmt=fmt, header=hdr, comments='')


def update_loggf(model, linelist, region='EWoptical'):
    '''Use correct loggf dependent on the model atmosphere and star

    Inputs
    ------
    model : str
      Model atmosphere to be used: "kurucz95" or "marcs"
    linelist : str
      Filename of the linelist to update
    region : str
      Spectral region: "EWoptical", "ABoptical", or "EWNIR"
    '''

    if model not in ['kurucz95', 'marcs']:
        raise IOError('Model not found: %s' % model)
    if region not in ['optical', 'NIR']:
        raise IOError('Region not found: %s' % region)

    df = pd.read_csv(linelist, skiprows=1, delimiter=r'\s+',
                     names=('wavelength', 'X', 'EP', 'loggf', 'EW'))

    if region == 'EWoptical':
        if model == 'kurucz95':
            loggf = 'Sousa2007_opt_kurucz.lst'
        elif model == 'marcs':
            loggf = 'Sousa2007_opt_marcs.lst'

    elif region == 'ABoptical':
        if model == 'kurucz95':
            loggf = 'Neves2009_opt_kurucz.lst'
        elif model == 'marcs':
            loggf = 'Neves2009_opt_marcs.lst'

    elif region == 'EWNIR':
        if model == 'kurucz95':
            loggf = 'Andreasen2016_nir_kurucz.lst'
        elif model == 'marcs':
            loggf = 'Andreasen2016_nir_marcs.lst'

    save_loggf(loggf, df, linelist)
