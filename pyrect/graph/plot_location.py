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
INPDIR = '/aos/home/kirk/projects/echo-class/echo/calibration/'

# Pickled height data
PRECIP = 'sgpxsaprppiI4.cloud.heights.pkl'
GROUND = 'sgpxsaprppiI4.ground.heights.pkl'
INSECTS = 'sgpxsaprppiI4.insect.heights.pkl'

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


def all_classes(filename, outdir=None, dpi=100, verbose=False):
    """
    """

    # Parse plot output directory
    if outdir is not None:
        filename = os.path.join(outdir, filename)

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

    fig = plt.figure(figsize=(8, 6))

    # (a) Heights
    axa = fig.add_subplot(111, xlim=(0, 16), ylim=(0, 1))
    axa.plot(data_p['bin centers'], data_p['normalized histogram'], 'k-',
             linewidth=2, label='Cloud')
    axa.plot(data_g['bin centers'], data_g['normalized histogram'], 'r-',
             linewidth=2, label='Ground')
    axa.plot(data_i['bin centers'], data_i['normalized histogram'], 'b-',
             linewidth=2, label='Insect')
    axa.xaxis.set_major_locator(MultipleLocator(2))
    axa.xaxis.set_minor_locator(MultipleLocator(1))
    axa.yaxis.set_major_locator(MultipleLocator(0.2))
    axa.yaxis.set_major_locator(MultipleLocator(0.1))
    axa.set_xlabel('(km)')
    axa.set_ylabel('Normalized Histogram')
    axa.set_title('Radar Gate Heights')
    axa.grid(which='major')

    # Add legend
    axa.legend(loc=[1.05, 0.4])

    # Save figure
    fig.savefig(filename, format='png', dpi=dpi, bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('filename', type=str, help=None)
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
        print 'filename = %s' % args.filename
        print 'outdir = %s' % args.outdir
        print 'dpi = %i' % args.dpi

    # Call desired plotting function
    all_classes(
        args.filename, outdir=args.outdir, dpi=args.dpi, verbose=args.verbose)
