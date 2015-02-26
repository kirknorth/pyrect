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


def histograms(pkl, outdir=None, dpi=100, verbose=False):
    """
    """

    # Parse plot output directory
    if outdir is None:
        outdir = ''

    # Read pickled data
    with open(pkl, 'rb') as fid:
        data = pickle.load(fid)

    if verbose:
        print [field for field in data.keys()]

    # Parse field textures
    phidp = data['Differential phase (PhiDP)']
    zdr = data['Differential reflectivity']
    rhohv = data['Cross correlation ratio (RHOHV)']
    refl = data['Reflectivity']

    fig = plt.figure(figsize=(18, 4))

    # (a) Differential phase texture
    axa = fig.add_subplot(141, xlim=(0, 180), ylim=(0, 1))
    axa.plot(phidp['bin centers'], phidp['normalized histogram'], 'k-',
             linewidth=2, label='Non-precip')
    axa.xaxis.set_major_locator(MultipleLocator(40))
    axa.xaxis.set_minor_locator(MultipleLocator(10))
    axa.yaxis.set_major_locator(MultipleLocator(0.2))
    axa.yaxis.set_major_locator(MultipleLocator(0.1))
    axa.set_xlabel('(deg)')
    axa.set_title('Differential phase texture')
    axa.grid(which='major')

    # (b) Differential reflectivity texture
    axb = fig.add_subplot(142, xlim=(0, 15), ylim=(0, 1))
    axb.plot(zdr['bin centers'], zdr['normalized histogram'], 'k-',
             linewidth=2, label='Non-precip')
    axb.xaxis.set_major_locator(MultipleLocator(3))
    axb.xaxis.set_minor_locator(MultipleLocator(1))
    axb.yaxis.set_major_locator(MultipleLocator(0.2))
    axb.yaxis.set_major_locator(MultipleLocator(0.1))
    axb.set_xlabel('(dBZ)')
    axb.set_title('Differential reflectivity texture')
    axb.grid(which='major')

    # (c) Copolar correlation texture
    axc = fig.add_subplot(143, xlim=(0, 0.6), ylim=(0, 1))
    axc.plot(rhohv['bin centers'], rhohv['normalized histogram'], 'k-',
             linewidth=2, label='Non-precip')
    axc.xaxis.set_major_locator(MultipleLocator(0.1))
    axc.xaxis.set_minor_locator(MultipleLocator(0.05))
    axc.yaxis.set_major_locator(MultipleLocator(0.2))
    axc.yaxis.set_major_locator(MultipleLocator(0.1))
    axc.set_title('Copolar correlation texture')
    axc.grid(which='major')

    # (d) Reflectivity texture
    axd = fig.add_subplot(144, xlim=(0, 30), ylim=(0, 1))
    axd.plot(refl['bin centers'], refl['normalized histogram'], 'k-',
             linewidth=2, label='Non-precip')
    axd.xaxis.set_major_locator(MultipleLocator(5))
    axd.xaxis.set_minor_locator(MultipleLocator(1))
    axd.yaxis.set_major_locator(MultipleLocator(0.2))
    axd.yaxis.set_major_locator(MultipleLocator(0.1))
    axd.set_xlabel('(dBZ)')
    axd.set_title('Reflectivity texture')
    axd.grid(which='major')

    # Add legend
    axd.legend(loc=[1.05, 0.4])

    # Save figure
    fname, fext = os.path.splitext(pkl)
    fname = '{}.png'.format(os.path.basename(fname))
    fig.savefig(os.path.join(outdir, fname), format='png', dpi=dpi,
                bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('pkl', type=str, help=None)
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
        print 'pkl = %s' % args.pkl
        print 'outdir = %s' % args.outdir
        print 'dpi = %i' % args.dpi

    # Call desired plotting function
    histograms(
        args.pkl, outdir=args.outdir, dpi=args.dpi, verbose=args.verbose)
