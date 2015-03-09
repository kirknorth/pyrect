#!/usr/bin/python

import os
import json
import getpass
import argparse
import numpy as np

from pyart.config import get_field_name

from clutter.texture import texture_fields

### GLOBAL VARIABLES ###
########################

# Define basic values and thresholds
MIN_NCP = 0.5
MIN_SAMPLE = 10
RAY_WINDOW = 5
GATE_WINDOW = 5
VCP_SWEEPS = 22
VCP_RAYS = 7920
MIN_SWEEP = 3
MAX_SWEEP = None

# Define bins and limits for texture histograms
BINS_PHIDP, LIMITS_PHIDP = 360, (0, 360)
BINS_ZDR, LIMITS_ZDR = 150, (0, 15)
BINS_RHOHV, LIMITS_RHOHV = 50, (0, 1)
BINS_REFL, LIMITS_REFL = 100, (0, 30)

# Define fields to exclude from radar object
EXCLUDE_FIELDS = [
    'radar_echo_classification',
    'spectrum_width',
    'corrected_reflectivity'
    ]

# Parse field names
ncp_field = get_field_name('normalized_coherent_power')
refl_field = get_field_name('reflectivity')
phidp_field = get_field_name('differential_phase')
zdr_field = get_field_name('differential_reflectivity')
rhohv_field = get_field_name('cross_correlation_ratio')


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=None)
    parser.add_argument('json', type=str, help=None)
    parser.add_argument('pickle', type=str, help=None)
    parser.add_argument('--inpdir', nargs='?', type=str, const='', default='',
                        help=None)
    parser.add_argument('--outdir', nargs='?', type=str, const='', default='',
                        help=None)
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
        print 'MIN_NCP = %.2f' % MIN_NCP
        print 'MIN_SAMPLE = %i' % MIN_SAMPLE
        print 'RAY_WINDOW = %i' % RAY_WINDOW
        print 'GATE_WINDOW = %i' % GATE_WINDOW
        print 'VCP_SWEEPS = %i' % VCP_SWEEPS
        print 'VCP_RAYS = %i' % VCP_RAYS
        print 'MIN_SWEEP = %s' % MIN_SWEEP
        print 'MAX_SWEEP = %s' % MAX_SWEEP

    # Compute histograms for specified fields
    phidp = texture_fields.histogram_from_json(
        args.json, phidp_field, inpdir=args.inpdir, ray_window=RAY_WINDOW,
        gate_window=GATE_WINDOW, min_sample=MIN_SAMPLE, bins=BINS_PHIDP,
        limits=LIMITS_PHIDP, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=ncp_field,
        verbose=args.verbose)

    zdr = texture_fields.histogram_from_json(
        args.json, zdr_field, inpdir=args.inpdir, ray_window=RAY_WINDOW,
        gate_window=GATE_WINDOW, min_sample=MIN_SAMPLE, bins=BINS_ZDR,
        limits=LIMITS_ZDR, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=ncp_field,
        verbose=args.verbose)

    rhohv = texture_fields.histogram_from_json(
        args.json, rhohv_field, inpdir=args.inpdir, ray_window=RAY_WINDOW,
        gate_window=GATE_WINDOW, min_sample=MIN_SAMPLE, bins=BINS_RHOHV,
        limits=LIMITS_RHOHV, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=ncp_field,
        verbose=args.verbose)

    refl = texture_fields.histogram_from_json(
        args.json, refl_field, inpdir=args.inpdir, ray_window=RAY_WINDOW,
        gate_window=GATE_WINDOW, min_sample=MIN_SAMPLE, bins=BINS_REFL,
        limits=LIMITS_REFL, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, min_sweep=MIN_SWEEP, max_sweep=MAX_SWEEP,
        exclude_fields=EXCLUDE_FIELDS, fill_value=None, ncp_field=ncp_field,
        verbose=args.verbose)

    # Pack histograms together
    histograms = [phidp, zdr, rhohv, refl]

    # Pickle texture histograms
    texture_fields._pickle_histograms(
        histograms, args.pickle, outdir=args.outdir)
