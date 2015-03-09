#!/usr/bin/python

import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MultipleLocator

from pyart.config import get_field_name


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


def histograms(precip, nonprecip, image, outdir=None, dpi=100, verbose=False):
    """
    """

    # Parse plot output directory
    if outdir is None:
        outdir = ''

    # Read pickled data
    with open(precip, 'rb') as fid:
        precip_data = pickle.load(fid)
    with open(nonprecip, 'rb') as fid:
        nonprecip_data = pickle.load(fid)

    if verbose:
        print 'Precipitation data: %s' % precip_data.keys()
        print 'Non-precipitation data: %s' % nonprecip_data.keys()

    # Parse moments data
    vdop_p = precip_data['velocity']
    vdop_np = nonprecip_data['velocity']
    sw_p = precip_data['spectrum_width']
    sw_np = nonprecip_data['spectrum_width']
    rhohv_p = precip_data['cross_correlation_ratio']
    rhohv_np = nonprecip_data['cross_correlation_ratio']

    fig = plt.figure(figsize=(18, 4))

    # (a) Doppler velocity
    axa = fig.add_subplot(141, xlim=(-20, 20), ylim=(0, 1))
    axa.plot(vdop_p['bin centers'], vdop_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axa.plot(vdop_np['bin centers'], vdop_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axa.xaxis.set_major_locator(MultipleLocator(5))
    axa.xaxis.set_minor_locator(MultipleLocator(1))
    axa.yaxis.set_major_locator(MultipleLocator(0.2))
    axa.yaxis.set_major_locator(MultipleLocator(0.1))
    axa.set_xlabel('(m/s)')
    axa.set_ylabel('Normalized Histogram')
    axa.set_title('Radial velocity')
    axa.grid(which='major')

    # (b) Spectrum width
    axb = fig.add_subplot(142, xlim=(0, 5), ylim=(0, 1))
    axb.plot(sw_p['bin centers'], sw_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axb.plot(sw_np['bin centers'], sw_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axb.xaxis.set_major_locator(MultipleLocator(1))
    axb.xaxis.set_minor_locator(MultipleLocator(0.5))
    axb.yaxis.set_major_locator(MultipleLocator(0.2))
    axb.yaxis.set_major_locator(MultipleLocator(0.1))
    axb.set_xlabel('(m/s)')
    axb.set_title('Spectrum width')
    axb.grid(which='major')

    # (c) Copolar correlation texture
    axc = fig.add_subplot(143, xlim=(0, 1), ylim=(0, 1))
    axc.plot(rhohv_p['bin centers'], rhohv_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axc.plot(rhohv_np['bin centers'], rhohv_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axc.xaxis.set_major_locator(MultipleLocator(0.2))
    axc.xaxis.set_minor_locator(MultipleLocator(0.1))
    axc.yaxis.set_major_locator(MultipleLocator(0.2))
    axc.yaxis.set_major_locator(MultipleLocator(0.1))
    axc.set_title('Copolar correlation coefficient')
    axc.grid(which='major')

    # Add legend
    axc.legend(loc=[1.03, 0.4])

    # Save figure
    fig.savefig(os.path.join(outdir, image), format='png', dpi=dpi,
                bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('precip', type=str, help=None)
    parser.add_argument('nonprecip', type=str, help=None)
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
        print 'precip = %s' % args.precip
        print 'nonprecip = %s' % args.nonprecip
        print 'image = %s' % args.image
        print 'outdir = %s' % args.outdir
        print 'dpi = %i' % args.dpi

    # Call desired plotting function
    histograms(
        args.precip, args.nonprecip, args.image, outdir=args.outdir,
        dpi=args.dpi, verbose=args.verbose)
