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
        filename, field, inpdir=None, bins=10, limits=None, min_ncp=0.3,
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


    return
