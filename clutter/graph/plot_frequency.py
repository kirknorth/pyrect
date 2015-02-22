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
NUM_SWEEPS = 22

# Define sweeps to be plotted
SWEEPS = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 18, 20, 21]

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

# Define color map
cmap = plt.get_cmap('jet')
norm = BoundaryNorm(np.arange(0, 1.1, 0.1), cmap.N)
ticks = np.arange(0, 1.1, 0.1)


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
    ax.set_title(title)

    return qm


def multipanel(filename, outdir, inpdir=None, dpi=50, verbose=False):
    """
    """

    # Parse input directory
    if inpdir is not None:
        filename = os.path.join(inpdir, filename)

    # Read pickled data
    with open(filename, 'rb') as fid:
        data = pickle.load(fid)

    # Parse radar data
    radar = data['last radar']

    if radar.nsweeps < NUM_SWEEPS:
        return

    if verbose:
        print 'Currently plotting file %s' % os.path.basename(filename)

    # Create figure instance
    subs = {'xlim': (-40, 40), 'ylim': (-40, 40)}
    figs = {'figsize': (52, 26)}
    fig, axes = plt.subplots(nrows=3, ncols=5, subplot_kw=subs, **figs)

    # Iterate over each sweep
    for k, ax in zip(SWEEPS, axes.flatten()):

        qm = _pcolormesh(radar, 'clutter_frequency_map', sweep=k, cmap=cmap,
                         norm=norm, ax=ax)

    # Format plot axes
    for ax in axes.flatten():
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(5))
        ax.yaxis.set_major_locator(MultipleLocator(10))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax.set_xlabel('Eastward Range from Radar (km)')
        ax.set_ylabel('Northward Range from Radar (km)')
        ax.grid(which='major')

    # Color bars
    cax = fig.add_axes([0.44, 0.05, 0.15, 0.01])
    cb = plt.colorbar(mappable=qm, cax=cax, extend='neither', ticks=ticks,
                      orientation='horizontal')
    cb.ax.set_ylabel('Prob.', rotation='horizontal', fontsize=24)
    cb.ax.yaxis.set_label_coords(1.07, 0.00)

    # Save figure
    fname, fext = os.path.splitext(filename)
    fname = '{}.png'.format(os.path.basename(fname))
    fig.savefig(os.path.join(outdir, fname), format='png', dpi=dpi,
                bbox_inches='tight')
    plt.close(fig)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('filename', type=str, help=None)
    parser.add_argument('outdir', type=str, help=None)
    parser.add_argument('-inpdir', nargs='?', type=str, const=None,
                        default=None, help=None)
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
        print 'inpdir = %s' % args.inpdir
        print 'dpi = %i' % args.dpi

    # Call desired plotting function
    multipanel(args.filename, args.outdir, inpdir=args.inpdir, dpi=args.dpi,
               verbose=args.verbose)



