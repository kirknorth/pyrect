#!/usr/bin/python

import os
import json
import getpass
import argparse
import numpy as np

from clutter.frequency import clutter_frequency

### GLOBAL VARIABLES ###
########################

# Define basic values and thresholds
MIN_NCP = 0.5
VCP_SWEEPS = 22
VCP_RAYS = 7920

# Define fields to exclude from radar object
EXCLUDE_FIELDS = [
    'reflectivity',
    'velocity',
    'cross_correlation_ratio',
    'radar_echo_classification',
    'corrected_reflectivity',
    'spectrum_width',
    'differential_phase',
    'differential_reflectivity'
    ]


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

    # Compute the clutter frequency map
    clutter = clutter_frequency.map_from_json(
        args.json, inpdir=args.inpdir, min_ncp=MIN_NCP, vcp_sweeps=VCP_SWEEPS,
        vcp_rays=VCP_RAYS, exclude_fields=EXCLUDE_FIELDS, ncp_field=None,
        verbose=args.verbose)

    # Pickle clutter frequency map
    clutter_frequency._pickle_map(clutter, args.pickle, outdir=args.outdir)
