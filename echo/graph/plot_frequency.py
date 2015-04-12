#!/usr/bin/python

import os
import pickle
import argparse
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MultipleLocator


### GLOBAL VARIABLES ###
########################

# Define the proper number of sweeps --> VCP to plot
VCP_SWEEPS = 22

# Define sweeps to be plotted
SWEEPS = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 21]

# Define color map
CMAP = plt.get_cmap('jet')
NORM = BoundaryNorm(np.arange(0, 1.1, 0.1), CMAP.N)
TICKS = np.arange(0, 1.1, 0.1)

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


def _pcolormesh(radar, field, sweep=0, cmap=None, norm=None, ax=None):
    """
    """

    if ax is None:
        ax = plt.gca()

    # Parse sweep data
    clutter = radar.get_field(sweep, field)
    azimuth = radar.get_azimuth(sweep)

    # Compute radar sweep coordinates
    AZI, RNG = np.meshgrid(
        np.radians(azimuth), radar.range['data'] / 1000.0, indexing='ij')
    X = RNG * np.sin(AZI)
    Y = RNG * np.cos(AZI)

    # Create plot
    qm = ax.pcolormesh(X, Y, clutter, cmap=cmap, norm=norm, shading='flat')

    # Create title
    title = '{} {:.1f} deg\n {}'.format(
       radar.metadata['instrument_name'], radar.fixed_angle['data'][sweep],
       radar.fields[field]['long_name'])
    ax.set_title(title, fontsize=18)

    return qm


def multipanel(pickle_file, outdir, inpdir=None, dpi=50, verbose=False):
    """
    """

    # Parse input directory
    if inpdir is not None:
        pickle_file = os.path.join(inpdir, pickle_file)

    # Read pickled data
    with open(pickle_file, 'rb') as fid:
        data = pickle.load(fid)

    # Parse radar data
    radar = data['last radar']

    if radar.nsweeps != VCP_SWEEPS:
        return

    if verbose:
        print 'Currently plotting file %s' % os.path.basename(pickle_file)

    # Create figure instance
    subs = {'xlim': (-40, 40), 'ylim': (-40, 40)}
    figs = {'figsize': (52, 26)}
    fig, ax = plt.subplots(nrows=3, ncols=5, subplot_kw=subs, **figs)
    ax = ax.flatten()

    # Iterate over each sweep
    for k, sweep in enumerate(SWEEPS):

        qm = _pcolormesh(
            radar, 'clutter_map', sweep=sweep, cmap=CMAP, norm=NORM, ax=ax[k])

    # Format plot axes
    for i in range(ax.size):
        ax[i].xaxis.set_major_locator(MultipleLocator(10))
        ax[i].xaxis.set_minor_locator(MultipleLocator(5))
        ax[i].yaxis.set_major_locator(MultipleLocator(10))
        ax[i].yaxis.set_minor_locator(MultipleLocator(5))
        ax[i].set_xlabel('Eastward Range from Radar (km)', fontsize=18)
        ax[i].set_ylabel('Northward Range from Radar (km)', fontsize=18)
        ax[i].grid(which='major')

    # Color bars
    cax = fig.add_axes([0.44, 0.05, 0.15, 0.01])
    cb = plt.colorbar(mappable=qm, cax=cax, extend='neither', ticks=TICKS,
                      orientation='horizontal')
    cb.ax.set_ylabel('Prob.', rotation='horizontal', fontsize=22)
    cb.ax.yaxis.set_label_coords(1.07, 0.00)

    # Save figure
    filename, ext = os.path.splitext(pickle_file)
    filename = '{}.png'.format(os.path.basename(filename))
    fig.savefig(os.path.join(outdir, filename), format='png', dpi=dpi,
                bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('pickle', type=str, help=None)
    parser.add_argument('outdir', type=str, help=None)
    parser.add_argument('--inpdir', nargs='?', type=str, const=None,
                        default=None, help=None)
    parser.add_argument('--dpi', nargs='?', type=int, const=50, default=50,
                        help=None)
    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True,
                        default=False, help=None)
    parser.add_argument('-db', '--debug', nargs='?', type=bool, const=True,
                        default=False, help=None)
    args = parser.parse_args()

    if args.debug:
        print 'pickle = {}'.format(args.pickle)
        print 'outdir = {}'.format(args.outdir)
        print 'inpdir = {}'.format(args.inpdir)
        print 'dpi = {}'.format(args.dpi)

    # Call desired plotting function
    multipanel(args.pickle, args.outdir, inpdir=args.inpdir, dpi=args.dpi,
               verbose=args.verbose)



