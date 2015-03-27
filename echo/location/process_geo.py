#!/usr/bin/python

import os
import json
import getpass
import argparse
import numpy as np

from pyart.io import read
from pyart.config import get_field_name

from echo.location import geo
from echo.correct import noise


### GLOBAL VARIABLES ###
########################

# Define basic values and thresholds
MIN_NCP = None
VCP_SWEEPS = 22
VCP_RAYS = None
MIN_SWEEP = None
MAX_SWEEP = None
MIN_RANGE = None
MAX_RANGE = None
REMOVE_SALT = True
SALT_WINDOW = (5, 5)
SALT_SAMPLE = 10

# Define bins and limits for location histograms
# Units should be in kilometers
BINS_HEIGHT, LIMITS_HEIGHT = 200, (0.0, 16.0)

# Define bins and limits for coherency (texuture) fields
BINS_VDOP_COHER, LIMITS_VDOP_COHER = 100, (0, 20)
BINS_SW_COHER, LIMITS_SW_COHER = 50, (0, 5)

# Define fields to exclude from radar object
EXCLUDE_FIELDS = [
    'reflectivity',
    'differential_phase',
    'differential_reflectivity',
    'cross_correlation_ratio',
    'total_power',
    'radar_echo_classification',
    'corrected_reflectivity',
    ]

# Parse field names
VDOP_FIELD = get_field_name('velocity')
SW_FIELD = get_field_name('spectrum_width')
NCP_FIELD = get_field_name('normalized_coherent_power')

# Create histogram dictionary
HIST_DICT = {
    'number of bins': BINS_HEIGHT,
    'limits': LIMITS_HEIGHT,
    'histogram counts': np.zeros(BINS_HEIGHT, dtype=np.float64),
    }


def _loop_over_dict(
        json_file, pickle_file, inpdir=None, outdir=None, verbose=False,
        debug=False):
    """
    """

    # Parse files from JSON
    with open(json_file, 'r') as fid:
        files = json.load(fid)

    if inpdir is not None:
        files = [os.path.join(inpdir, f) for f in files]

    # Loop over all files
    for f in files:

        # Read radar data
        radar = read(f, exclude_fields=EXCLUDE_FIELDS)

        if VCP_SWEEPS is not None and radar.nsweeps != VCP_SWEEPS:
            continue
        if VCP_RAYS is not None and radar.nrays != VCP_RAYS:
            continue

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Determine significant detection of the radar
        gatefilter = noise.velocity_coherency(
            radar, gatefilter=None, num_bins=BINS_VDOP_COHER,
            limits=LIMITS_VDOP_COHER, texture_window=(3, 3),
            texture_sample=5, min_sigma=None, max_sigma=None, nyquist=None,
            rays_wrap_around=False, remove_salt=REMOVE_SALT,
            salt_window=SALT_WINDOW, salt_sample=SALT_SAMPLE, fill_value=None,
            verbose=verbose)
        gatefilter = noise.spectrum_width_coherency(
            radar, gatefilter=gatefilter, num_bins=BINS_SW_COHER,
            limits=LIMITS_SW_COHER, texture_window=(3, 3), texture_sample=5,
            min_sigma=None, max_sigma=None, rays_wrap_around=False,
            remove_salt=REMOVE_SALT, salt_window=SALT_WINDOW,
            salt_sample=SALT_SAMPLE, fill_value=None, verbose=True)
        gatefilter = noise.significant_detection(
            radar, gatefilter=gatefilter, remove_salt=REMOVE_SALT,
            salt_window=SALT_WINDOW, salt_sample=SALT_SAMPLE, min_ncp=MIN_NCP,
            detect_field=None, verbose=verbose)

        # Compute histogram counts for each field
        geo.height_histogram_from_radar(
            radar, HIST_DICT, gatefilter=gatefilter, min_ncp=MIN_NCP,
            min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP, min_range=MIN_RANGE,
            max_range=MAX_RANGE, fill_value=None, ncp_field=NCP_FIELD,
            verbose=verbose, debug=debug)

    # Parse bin edges and histogram counts
    bin_edges = HIST_DICT['bin edges']
    counts = HIST_DICT['histogram counts']

    # Compute normalized histogram and probability density
    # Add these to the histogram dictionary
    counts_norm = counts.astype(np.float64) / counts.max()
    pdf = counts_norm / np.sum(counts_norm * np.diff(bin_edges))
    HIST_DICT['normalized histogram'] = counts_norm
    HIST_DICT['probability density'] = pdf

    # Include other parameters in the histogram dictionary
    HIST_DICT['radar files'] = files
    HIST_DICT['min sweep'] = MIN_SWEEP
    HIST_DICT['max sweep'] = MAX_SWEEP
    HIST_DICT['min range'] = MIN_RANGE
    HIST_DICT['max range'] = MAX_RANGE
    HIST_DICT['sweeps in VCP'] = VCP_SWEEPS
    HIST_DICT['rays in VCP'] = VCP_RAYS
    HIST_DICT['min NCP'] = MIN_NCP

    # Pickle histogram data
    geo._pickle_histograms(
        HIST_DICT, pickle_file, outdir=outdir)

    return


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('json', type=str, help=None)
    parser.add_argument('pickle', type=str, help=None)
    parser.add_argument('--inpdir', nargs='?', type=str, const=None,
                        default=None, help=None)
    parser.add_argument('--outdir', nargs='?', type=str, const=None,
                        default=None, help=None)
    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True,
                        default=False, help=None)
    parser.add_argument('-db', '--debug', nargs='?', type=bool, const=True,
                        default=False, help=None)
    args = parser.parse_args()

    if args.debug:
        print 'json = %s' % args.json
        print 'pickle = %s' % args.pickle
        print 'inpdir = %s' % args.inpdir
        print 'outdir = %s' % args.outdir

    if args.verbose:
        print 'MIN_NCP = {}'.format(MIN_NCP)
        print 'VCP_SWEEPS = {}'.format(VCP_SWEEPS)
        print 'VCP_RAYS = {}'.format(VCP_RAYS)
        print 'MIN_SWEEP = {}'.format(MIN_SWEEP)
        print 'MAX_SWEEP = {}'.format(MAX_SWEEP)
        print 'MIN_RANGE = {} km'.format(MIN_RANGE)
        print 'MAX_RANGE = {} km'.format(MAX_RANGE)
        print 'REMOVE_SALT = {}'.format(REMOVE_SALT)
        print 'SALT_WINDOW = {}'.format(SALT_WINDOW)
        print 'SALT_SAMPLE = {}'.format(SALT_SAMPLE)

    # Use desired processing method
    _loop_over_dict(
        args.json, args.pickle, inpdir=args.inpdir, outdir=args.outdir,
        verbose=args.verbose, debug=args.debug)

