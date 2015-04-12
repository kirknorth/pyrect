#!/usr/bin/python

import os
import json
import getpass
import argparse
import numpy as np

from pyart.io import read
from pyart.config import get_field_name

from echo.moment import moment_fields
from echo.correct import noise


### GLOBAL VARIABLES ###
########################

# Define basic values and thresholds
MIN_NCP = None
VCP_SWEEPS = 22
VCP_RAYS = None
MIN_SWEEP = 2
MAX_SWEEP = None
MIN_RANGE = 3.0
MAX_RANGE = None
REMOVE_SALT = True
SALT_WINDOW = (5, 5)
SALT_SAMPLE = 10

# Define bins and limits for moment histograms
BINS_REFL, LIMITS_REFL = 200, (-40.0, 50.0)
BINS_VDOP, LIMITS_VDOP = 200, (-20.0, 20.0)
BINS_SW, LIMITS_SW = 100, (0, 10)
BINS_RHOHV, LIMITS_RHOHV = 50, (0, 1)
BINS_ZDR, LIMITS_ZDR = 250, (-20, 30)
BINS_NCP, LIMITS_NCP = 100, (0, 1)

# Define bins and limits for coherency (texuture) fields
BINS_VDOP_COHER, LIMITS_VDOP_COHER = 100, (0, 20)
BINS_SW_COHER, LIMITS_SW_COHER = 50, (0, 5)

# Define fields to exclude from radar object
EXCLUDE_FIELDS = [
    'radar_echo_classification',
    'corrected_reflectivity',
    'differential_phase',
    ]

# Parse field names
REFL_FIELD = get_field_name('reflectivity')
VDOP_FIELD = get_field_name('velocity')
SW_FIELD = get_field_name('spectrum_width')
RHOHV_FIELD = get_field_name('cross_correlation_ratio')
ZDR_FIELD = get_field_name('differential_reflectivity')
NCP_FIELD = get_field_name('normalized_coherent_power')
PHIDP_FIELD = get_field_name('differential_phase')

# Create histogram dictionary
HIST_DICT = {
    REFL_FIELD: {'number of bins': BINS_REFL,
                 'limits': LIMITS_REFL,
                 'histogram counts': np.zeros(BINS_REFL, dtype=np.float64),
                 },
    VDOP_FIELD: {'number of bins': BINS_VDOP,
                 'limits': LIMITS_VDOP,
                 'histogram counts': np.zeros(BINS_VDOP, dtype=np.float64),
                 },
    SW_FIELD: {'number of bins': BINS_SW,
               'limits': LIMITS_SW,
               'histogram counts': np.zeros(BINS_SW, dtype=np.float64),
               },
    RHOHV_FIELD: {'number of bins': BINS_RHOHV,
                  'limits': LIMITS_RHOHV,
                  'histogram counts': np.zeros(BINS_RHOHV, dtype=np.float64),
                  },
    ZDR_FIELD: {'number of bins': BINS_ZDR,
                'limits': LIMITS_ZDR,
                'histogram counts': np.zeros(BINS_ZDR, dtype=np.float64),
                },
    NCP_FIELD: {'number of bins': BINS_NCP,
                'limits': LIMITS_NCP,
                'histogram counts': np.zeros(BINS_NCP, dtype=np.float64),
                },
    }


def _loop_over_fields(json, pickle, inpdir=None, outdir=None, verbose=False):
    """
    """

    # Reflectivity
    if verbose:
        print 'Processing reflectivity'
    refl = moment_fields.histogram_from_json(
        json, REFL_FIELD, inpdir=inpdir, bins=BINS_REFL,
        limits=LIMITS_REFL, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Doppler velocity
    if verbose:
        print 'Processing Doppler velocity'
    vdop = moment_fields.histogram_from_json(
        json, VDOP_FIELD, inpdir=inpdir, bins=BINS_VDOP,
        limits=LIMITS_VDOP, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Spectrum width
    if verbose:
        print 'Processing spectrum width'
    sw = moment_fields.histogram_from_json(
        json, SW_FIELD, inpdir=inpdir, bins=BINS_SW,
        limits=LIMITS_SW, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Copolar correlation
    if verbose:
        print 'Processing copolar correlation'
    rhohv = moment_fields.histogram_from_json(
        json, RHOHV_FIELD, inpdir=inpdir, bins=BINS_RHOHV,
        limits=LIMITS_RHOHV, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Differential reflectivity
    if verbose:
        print 'Processing differential reflectivity'
    zdr = moment_fields.histogram_from_json(
        json, ZDR_FIELD, inpdir=inpdir, bins=BINS_ZDR,
        limits=LIMITS_ZDR, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Coherent power
    if verbose:
        print 'Processing normalized coherent power'
    ncp = moment_fields.histogram_from_json(
        json, NCP_FIELD, inpdir=inpdir, bins=BINS_NCP,
        limits=LIMITS_NCP, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=NCP_FIELD,
        verbose=verbose)

    # Pack histograms together
    histograms = [refl, vdop, sw, rhohv, zdr, ncp]

    # Pickle moment histograms
    moment_fields._pickle_histograms(
        histograms, pickle, outdir=outdir)

    return


def _loop_over_dict(json_file, pickle_file, inpdir=None, outdir=None,
                    verbose=False):
    """
    """

    # Parse files from JSON file
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
            verbose=True)
        gatefilter = noise.spectrum_width_coherency(
            radar, gatefilter=gatefilter, num_bins=BINS_SW_COHER,
            limits=LIMITS_SW_COHER, texture_window=(3, 3), texture_sample=5,
            min_sigma=None, max_sigma=None, rays_wrap_around=False,
            remove_salt=REMOVE_SALT, salt_window=SALT_WINDOW,
            salt_sample=SALT_SAMPLE, fill_value=None, verbose=True)
        gatefilter = noise.significant_detection(
            radar, gatefilter=gatefilter, remove_salt=REMOVE_SALT,
            salt_window=SALT_WINDOW, salt_sample=SALT_SAMPLE, min_ncp=MIN_NCP,
            detect_field=None, verbose=True)

        # Compute histogram counts for each field
        moment_fields.histograms_from_radar(
            radar, HIST_DICT, gatefilter=gatefilter, min_ncp=MIN_NCP,
            min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP, min_range=MIN_RANGE,
            max_range=MAX_RANGE, fill_value=None, ncp_field=NCP_FIELD,
            verbose=verbose)

    # Normalize histograms for each field and compute probability densities
    for field in HIST_DICT:

        # Parse bin edges and histogram counts
        bin_edges = HIST_DICT[field]['bin edges']
        counts = HIST_DICT[field]['histogram counts']

        # Compute normalized histogram and probability density
        # Add these to the histogram dictionary
        counts_norm = counts.astype(np.float64) / counts.max()
        pdf = counts_norm / np.sum(counts_norm * np.diff(bin_edges))
        HIST_DICT[field]['normalized histogram'] = counts_norm
        HIST_DICT[field]['probability density'] = pdf

        # Include other parameters in the histogram dictionary
        HIST_DICT[field]['radar files'] = files
        HIST_DICT[field]['min sweep'] = MIN_SWEEP
        HIST_DICT[field]['max sweep'] = MAX_SWEEP
        HIST_DICT[field]['min range'] = MIN_RANGE
        HIST_DICT[field]['max range'] = MAX_RANGE
        HIST_DICT[field]['sweeps in VCP'] = VCP_SWEEPS
        HIST_DICT[field]['rays in VCP'] = VCP_RAYS
        HIST_DICT[field]['min NCP'] = MIN_NCP

    # Pickle histogram data
    moment_fields._pickle_histograms(
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
    parser.add_argument('--method', nargs='?', type=str, const='dict',
                        default='dict', help=None)
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
        print 'method = %s' % args.method

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
    if args.method.upper() == 'DICT':
        _loop_over_dict(
            args.json, args.pickle, inpdir=args.inpdir, outdir=args.outdir,
            verbose=args.verbose)

    elif args.method.upper() == 'FIELDS':
        _loop_over_fields(
            args.json, args.pickle, inpdir=args.inpdir, outdir=args.outdir,
            verbose=args.verbose)

    else:
        raise ValueError('Unsupported processing method')

