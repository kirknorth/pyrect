#!/usr/bin/python

import os
import json
import getpass
import argparse
import numpy as np

from echo.frequency import nonprecip


### GLOBAL VARIABLES ###
########################

# Define basic values and thresholds
# Define basic values and thresholds
MIN_NCP = 0.3
VCP_SWEEPS = 22
VCP_RAYS = None
VCP_GATES = None
MIN_SWEEP = None
MAX_SWEEP = None
MIN_RANGE = None
MAX_RANGE = None
USE_FILTER = False
TEXTURE_WINDOW = (3, 3)
TEXTURE_SAMPLE = 5
REMOVE_SALT = True
SALT_WINDOW = (5, 5)
SALT_SAMPLE = 10

# Define bins and limits for coherency fields
BINS_VDOP, LIMITS_VDOP = 100, (0, 20)
BINS_SW, LIMITS_SW = 50, (0, 5)

# Define fields to exclude from radar object
EXCLUDE_FIELDS = [
    'reflectivity',
    'total_power',
    'cross_correlation_ratio',
    'radar_echo_classification',
    'corrected_reflectivity',
    'differential_phase',
    'differential_reflectivity',
    ]


def _map_from_json(
        json_file, pickle_file, inpdir=None, outdir=None, verbose=False,
        debug=False):
    """
    """

    # Compute non-precipitating frequency map from JSON file
    map_dict = nonprecip.map_from_json(
        json_file, inpdir=inpdir, vcp_sweeps=VCP_SWEEPS, vcp_rays=VCP_RAYS,
        vcp_gates=VCP_GATES, min_ncp=MIN_NCP, use_filter=USE_FILTER,
        texture_window=TEXTURE_WINDOW, texture_sample=TEXTURE_SAMPLE,
        vdop_bins=BINS_VDOP, vdop_limits=LIMITS_VDOP, sw_bins=BINS_SW,
        sw_limits=LIMITS_SW, remove_salt=REMOVE_SALT, salt_window=SALT_WINDOW,
        salt_sample=SALT_SAMPLE, exclude_fields=EXCLUDE_FIELDS, ncp_field=None,
        debug=debug, verbose=verbose)

    # Pickle non-precipitating frequency map
    nonprecip._pickle_map(map_dict, pickle_file, outdir=outdir)

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
    parser.add_argument('--method', nargs='?', type=str, const='json',
                        default='json', help=None)
    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True,
                        default=False, help=None)
    parser.add_argument('-db', '--debug', nargs='?', type=bool, const=True,
                        default=False, help=None)
    args = parser.parse_args()

    if args.debug:
        print 'json = {}'.format(args.json)
        print 'pickle = {}'.format(args.pickle)
        print 'inpdir = {}'.format(args.inpdir)
        print 'outdir = {}'.format(args.outdir)
        print 'method = {}'.format(args.method)

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

    # Use desired processing function
    if args.method.upper() == 'JSON':
        _map_from_json(
            args.json, args.pickle, inpdir=args.inpdir, outdir=args.outdir,
            verbose=args.verbose)

    else:
        raise ValueError('Unsupported processing function')
