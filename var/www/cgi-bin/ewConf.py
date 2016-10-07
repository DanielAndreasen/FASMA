#!/home/daniel/Software/anaconda3/bin/python

# Import modules for CGI handling
import cgi, cgitb
from aresDriver import aresdriver
import pandas as pd

def cgi2dict(form):
    """Convert the form from cgi.FieldStorage to a python dictionary"""
    params = {}
    for key in form.keys():
        params[key] = form[key].value
    return params


def ares(form):
    """Create the configuration file for running the ARES driver"""
    def rejt_from_snr(snr):
        """Calculate rejt from SNR"""
        return 1-1/snr

    if form['continuum'] == 'snr':
        rejt = rejt_from_snr(form['snr'])
    elif form['continuum'] == 'rejt':
        rejt = form['rejt']
        rejt = 0.99999 if rejt > 0.99999 else rejt
    rvmask = False if form['rv'] == 0 else form['rv']

    # Make the StarMe_ares.cfg file from Gooey
    if form['linelist'] == 'Optical (parameters)':
        fout = 'Sousa2007_opt.lst'
    elif form['linelist'] == 'NIR (parameters)':
        fout = 'Andreasen2016_nir.lst'
    elif form['linelist'] == 'Optical (abundances)':
        fout = 'Neves2009_opt.lst'
    fout += ' %s' % form['spectrum']
    fout += ' lambdai:%s,lambdaf:%s,smoothder:%s' % (form['w0'], form['wf'], form['smooth'])
    fout += ',space:%s,lineresol:%s' % (form['space'], form['lineresol'])
    fout += ',miniline:%s,EWcut:%s' % (form['miniline'], form['EWcut'])
    if form['continuum'] in ['snr', 'rejt']:
        fout += ',rejt:%s' % rejt
    if rvmask:
        fout += ',rvmask:%s' % rvmask
    if 'force' in form.keys():
        fout += ',force'

    with open('/tmp/StarMe_ares.cfg', 'w') as f:
        f.writelines(fout)

    aresdriver('/tmp/StarMe_ares.cfg')
    linelist_out = 'linelist/%s.moog' % form['spectrum'].rpartition('.')[0]


if __name__ == '__main__':
    cgitb.enable()
    # Create instance of FieldStorage
    form = cgi2dict(cgi.FieldStorage())
    ares(form)
    print "Content-type: text/html\n\n"
    print "<html>"
    print "<head>"
    print "<title>FASMA - ARES driver</title>"
    print "</head>"
    print "<body>"
    for key in form.iterkeys():
        print "<p>%s: %s</p>" % (key, form[key])
    print "<p>%s</p>" % form
    print "<p>Pandas version: %s</p>" % pd.__version__
    print "</body>"
    print "</html>"
