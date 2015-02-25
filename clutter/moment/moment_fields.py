"""
clutter.moment.moment_fields
============================

"""

import os
import pickle
import numpy as np


def _pickle_histograms(histograms, filename, outdir=None):
    """
    """


def histogram_from_json(
        filename, field, inpdir=None, bins=10, limits=None, min_ncp=0.5,
        vcp_sweeps=22, min_sweep=None, max_sweep=None, exclude_fields=None,
        fill_value=None, ncp_field=None):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Parse files from JSON file
    with open(filename, 'r') as fid:
        files = json.load(fid)

    # Append input directory if given
    if inpdir is not None:
        files = [os.path.join(inpdir, f) for f in files]

    if verbose:
        print 'Total number of radar files to process = %i' % len(files)

    # Loop over all files
    histogram = np.zeros(bins, dtype=np.float64)
    for f in files:

        # Read radar data
        radar = read(f, exclude_fields=exclude_fields)

        if radar.nsweeps != vcp_sweeps:
            continue

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Parse data and compute histogram
        data = radar.fields[field]['data']
        hist, edges = np.histogram(
            data.compressed(), bins=bins, range=limits, normed=False,
            weights=None, density=False)
        histogram += hist

    return {
        'field': radar.fields[field]['long_name'],
        'histogram': histogram,
        'normalized histogram': histogram / histogram.max(),
        'number of bins': bins,
        'limits': limits,
        'bin edges': edges,
        'bin centers': edges[:-1] + np.diff(edges) / 2.0,
        'radar files': [os.path.basename(f) for f in files],
        'min sweep': min_sweep,
        'max sweep': max_sweep,
        'min normalized coherent power': min_ncp,
        'sweeps in VCP': vcp_sweeps,
        }
