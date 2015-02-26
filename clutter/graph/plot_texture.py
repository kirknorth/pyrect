#!/usr/bin/python

import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MultipleLocator


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

    # Parse field textures
    phidp_p = precip_data['Differential phase (PhiDP)']
    phidp_np = nonprecip_data['Differential phase (PhiDP)']
    zdr_p = precip_data['Differential reflectivity']
    zdr_np = nonprecip_data['Differential reflectivity']
    rhohv_p = precip_data['Cross correlation ratio (RHOHV)']
    rhohv_np = nonprecip_data['Cross correlation ratio (RHOHV)']
    refl_p = precip_data['Reflectivity']
    refl_np = nonprecip_data['Reflectivity']

    fig = plt.figure(figsize=(18, 4))

    # (a) Differential phase texture
    axa = fig.add_subplot(141, xlim=(0, 180), ylim=(0, 1))
    axa.plot(phidp_p['bin centers'], phidp_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axa.plot(phidp_np['bin centers'], phidp_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axa.xaxis.set_major_locator(MultipleLocator(40))
    axa.xaxis.set_minor_locator(MultipleLocator(10))
    axa.yaxis.set_major_locator(MultipleLocator(0.2))
    axa.yaxis.set_major_locator(MultipleLocator(0.1))
    axa.set_xlabel('(deg)')
    axa.set_ylabel('Normalized Histogram')
    axa.set_title('Differential phase texture')
    axa.grid(which='major')

    # (b) Differential reflectivity texture
    axb = fig.add_subplot(142, xlim=(0, 12), ylim=(0, 1))
    axb.plot(zdr_p['bin centers'], zdr_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axb.plot(zdr_np['bin centers'], zdr_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axb.xaxis.set_major_locator(MultipleLocator(2))
    axb.xaxis.set_minor_locator(MultipleLocator(1))
    axb.yaxis.set_major_locator(MultipleLocator(0.2))
    axb.yaxis.set_major_locator(MultipleLocator(0.1))
    axb.set_xlabel('(dBZ)')
    axb.set_title('Differential reflectivity texture')
    axb.grid(which='major')

    # (c) Copolar correlation texture
    axc = fig.add_subplot(143, xlim=(0, 0.5), ylim=(0, 1))
    axc.plot(rhohv_p['bin centers'], rhohv_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axc.plot(rhohv_np['bin centers'], rhohv_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axc.xaxis.set_major_locator(MultipleLocator(0.1))
    axc.xaxis.set_minor_locator(MultipleLocator(0.05))
    axc.yaxis.set_major_locator(MultipleLocator(0.2))
    axc.yaxis.set_major_locator(MultipleLocator(0.1))
    axc.set_title('Copolar correlation texture')
    axc.grid(which='major')

    # (d) Reflectivity texture
    axd = fig.add_subplot(144, xlim=(0, 25), ylim=(0, 1))
    axd.plot(refl_p['bin centers'], refl_p['normalized histogram'], 'k-',
             linewidth=2, label='Precip')
    axd.plot(refl_np['bin centers'], refl_np['normalized histogram'], 'r-',
             linewidth=2, label='Non-precip')
    axd.xaxis.set_major_locator(MultipleLocator(5))
    axd.xaxis.set_minor_locator(MultipleLocator(1))
    axd.yaxis.set_major_locator(MultipleLocator(0.2))
    axd.yaxis.set_major_locator(MultipleLocator(0.1))
    axd.set_xlabel('(dBZ)')
    axd.set_title('Reflectivity texture')
    axd.grid(which='major')

    # Add legend
    axd.legend(loc=[1.03, 0.4])

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
