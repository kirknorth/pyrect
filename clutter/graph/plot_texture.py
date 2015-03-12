#!/usr/bin/python

import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator

from pyart.config import get_field_name


### GLOBAL VARIABLES ###

# Input directory to moments data
INPDIR = '/aos/home/kirk/projects/clutter-classification/clutter/calibration'

# Pickled moments data
PRECIP = 'sgpxsaprppiI4.precip.textures.pkl'
GROUND = 'sgpxsaprppiI4.ground.textures.pkl'
INSECTS = 'sgpxsaprppiI4.insects.textures.pkl'

# Parse field names
refl_field = get_field_name('reflectivity')
vdop_field = get_field_name('velocity')
sw_field = get_field_name('spectrum_width')
rhv_field = get_field_name('cross_correlation_ratio')
zdr_field = get_field_name('differential_reflectivity')
phi_field = get_field_name('differential_phase')
ncp_field = get_field_name('normalized_coherent_power')

### Set figure parameters ###
rcParams['axes.linewidth'] = 1.5
rcParams['xtick.major.size'] = 4
rcParams['xtick.major.width'] = 1
rcParams['xtick.minor.size'] = 2
rcParams['xtick.minor.width'] = 1
rcParams['ytick.major.size'] = 4
rcParams['ytick.major.width'] = 1
rcParams['ytick.minor.size'] = 2
rcParams['ytick.minor.width'] = 1


def all_classes(image, outdir=None, dpi=100, verbose=False):
    """
    """

    # Parse plot output directory
    if outdir is None:
        outdir = ''

    # Read pickled data
    with open(os.path.join(INPDIR, PRECIP), 'rb') as fid:
        data_p = pickle.load(fid)
    with open(os.path.join(INPDIR, GROUND), 'rb') as fid:
        data_g = pickle.load(fid)
    with open(os.path.join(INPDIR, INSECTS), 'rb') as fid:
        data_i = pickle.load(fid)

    if verbose:
        print 'Precipitation data: %s' % data_p.keys()
        print 'Ground clutter data: %s' % data_g.keys()
        print 'Insects data: %s' % data_i.keys()

    # Parse moments data
    refl_p = data_p[refl_field]
    refl_g = data_g[refl_field]
    refl_i = data_i[refl_field]
    vdop_p = data_p[vdop_field]
    vdop_g = data_g[vdop_field]
    vdop_i = data_i[vdop_field]
    sw_p = data_p[sw_field]
    sw_g = data_g[sw_field]
    sw_i = data_i[sw_field]
    rhohv_p = data_p[rhv_field]
    rhohv_g = data_g[rhv_field]
    rhohv_i = data_i[rhv_field]
    zdr_p = data_p[zdr_field]
    zdr_g = data_g[zdr_field]
    zdr_i = data_i[zdr_field]
    phi_p = data_p[phi_field]
    phi_g = data_g[phi_field]
    phi_i = data_i[phi_field]
    ncp_p = data_p[ncp_field]
    ncp_g = data_g[ncp_field]
    ncp_i = data_i[ncp_field]

    fig = plt.figure(figsize=(20, 6))

    # (a) Differential phase texture
    axa = fig.add_subplot(241, xlim=(0, 180), ylim=(0, 1))
    axa.plot(phi_p['bin centers'], phi_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axa.plot(phi_g['bin centers'], phi_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axa.plot(phi_i['bin centers'], phi_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axa.xaxis.set_major_locator(MultipleLocator(40))
    axa.xaxis.set_minor_locator(MultipleLocator(10))
    axa.yaxis.set_major_locator(MultipleLocator(0.2))
    axa.yaxis.set_major_locator(MultipleLocator(0.1))
    axa.set_xlabel('(deg)')
    axa.set_ylabel('Normalized Histogram')
    axa.set_title('Differential phase texture')
    axa.grid(which='major')

    # (b) Differential reflectivity texture
    axb = fig.add_subplot(242, xlim=(0, 12), ylim=(0, 1))
    axb.plot(zdr_p['bin centers'], zdr_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axb.plot(zdr_g['bin centers'], zdr_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axb.plot(zdr_i['bin centers'], zdr_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axb.xaxis.set_major_locator(MultipleLocator(2))
    axb.xaxis.set_minor_locator(MultipleLocator(1))
    axb.yaxis.set_major_locator(MultipleLocator(0.2))
    axb.yaxis.set_major_locator(MultipleLocator(0.1))
    axb.set_xlabel('(dBZ)')
    axb.set_title('Differential reflectivity texture')
    axb.grid(which='major')

    # (c) Copolar correlation texture
    axc = fig.add_subplot(243, xlim=(0, 0.5), ylim=(0, 1))
    axc.plot(rhohv_p['bin centers'], rhohv_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axc.plot(rhohv_g['bin centers'], rhohv_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axc.plot(rhohv_i['bin centers'], rhohv_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axc.xaxis.set_major_locator(MultipleLocator(0.1))
    axc.xaxis.set_minor_locator(MultipleLocator(0.05))
    axc.yaxis.set_major_locator(MultipleLocator(0.2))
    axc.yaxis.set_major_locator(MultipleLocator(0.1))
    axc.set_title('Copolar correlation texture')
    axc.grid(which='major')

    # (d) Reflectivity texture
    axd = fig.add_subplot(244, xlim=(0, 25), ylim=(0, 1))
    axd.plot(refl_p['bin centers'], refl_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axd.plot(refl_g['bin centers'], refl_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axd.plot(refl_i['bin centers'], refl_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axd.xaxis.set_major_locator(MultipleLocator(5))
    axd.xaxis.set_minor_locator(MultipleLocator(1))
    axd.yaxis.set_major_locator(MultipleLocator(0.2))
    axd.yaxis.set_major_locator(MultipleLocator(0.1))
    axd.set_xlabel('(dBZ)')
    axd.set_title('Reflectivity texture')
    axd.grid(which='major')

    # (e) Doppler velocity texture
    axe = fig.add_subplot(245, xlim=(0, 15), ylim=(0, 1))
    axe.plot(vdop_p['bin centers'], vdop_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axe.plot(vdop_g['bin centers'], vdop_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axe.plot(vdop_i['bin centers'], vdop_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axe.xaxis.set_major_locator(MultipleLocator(3))
    axe.xaxis.set_minor_locator(MultipleLocator(1))
    axe.yaxis.set_major_locator(MultipleLocator(0.2))
    axe.yaxis.set_major_locator(MultipleLocator(0.1))
    axe.set_xlabel('(m/s)')
    axe.set_title('Radial velocity texture')
    axe.grid(which='major')

    # (f) Spectrum width texture
    axf = fig.add_subplot(246, xlim=(0, 10), ylim=(0, 1))
    axf.plot(sw_p['bin centers'], sw_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axf.plot(sw_g['bin centers'], sw_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axf.plot(sw_i['bin centers'], sw_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axf.xaxis.set_major_locator(MultipleLocator(2))
    axf.xaxis.set_minor_locator(MultipleLocator(1))
    axf.yaxis.set_major_locator(MultipleLocator(0.2))
    axf.yaxis.set_major_locator(MultipleLocator(0.1))
    axf.set_xlabel('(m/s)')
    axf.set_title('Spectrum width texture')
    axf.grid(which='major')

    # (g) Normalized coherent power texture
    axg = fig.add_subplot(247, xlim=(0, 1), ylim=(0, 1))
    axg.plot(ncp_p['bin centers'], ncp_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axg.plot(ncp_g['bin centers'], ncp_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axg.plot(ncp_i['bin centers'], ncp_i['normalized histogram'], 'b-',
             linewidth=2, label='Insects')
    axg.xaxis.set_major_locator(MultipleLocator(2))
    axg.xaxis.set_minor_locator(MultipleLocator(1))
    axg.yaxis.set_major_locator(MultipleLocator(0.2))
    axg.yaxis.set_major_locator(MultipleLocator(0.1))
    axg.set_xlabel('')
    axg.set_title('Normalized coherent power texture')
    axg.grid(which='major')

    # Add legend
    axd.legend(loc=[1.05, -0.10])

    # Save figure
    fig.savefig(os.path.join(outdir, image), format='png', dpi=dpi,
                bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('image', type=str, help=None)
    parser.add_argument('--outdir', nargs='?', type=str, const='', default='',
                        help=None)
    parser.add_argument('--dpi', nargs='?', type=int, const=50, default=50,
                        help=None)
    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True,
                        default=False, help=None)
    parser.add_argument('-db', '--debug', nargs='?', type=bool, const=True,
                        default=False, help=None)
    args = parser.parse_args()

    if args.debug:
        print 'image = %s' % args.image
        print 'outdir = %s' % args.outdir
        print 'dpi = %i' % args.dpi

    # Call desired plotting function
    histograms(
        args.image, outdir=args.outdir, dpi=args.dpi, verbose=args.verbose)
