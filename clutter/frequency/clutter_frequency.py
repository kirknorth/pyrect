"""
clutter.frequency.clutter_frequency
===================================

"""

import os
import re
import json
import pickle
import numpy as np

from datetime import datetime

from pyart.io import read
from pyart.config import get_field_name


def _pickle_map(clutter_map, filename, outdir=None):
    """
    Pickle a clutter frequency (probability) map.
    """

    # Parse output directory
    if outdir is not None:
        filename = os.path.join(outdir, filename)

    with open(os.path.join(outdir, filename), 'wb') as fid:
        pickle.dump(clutter_map, fid, protocol=pickle.HIGHEST_PROTOCOL)

    return


def map_date_range(start, stop, inpdir, stamp, date_str='[0-9]{12}',
                   date_fmt='%y%m%d%H%M%S', min_ncp=0.5, vcp_sweeps=None,
                   vcp_rays=None, exclude_fields=None, ncp_field=None,
                   debug=False, verbose=False):
    """
    Compute the clutter frequency (probability) map within the specified date
    range. The start and stop times should define a non-precipitating time
    period where (most) echoes present must be, by definition, clutter.

    Parameters
    ----------

    Optional Parameters
    -------------------

    Returns
    -------
    """

    # Parse field names
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Get all files with stamp in directory
    files = [os.path.join(inpdir, f) for f in sorted(os.listdir(inpdir))
             if stamp in f]

    if verbose:
        print 'Total number of radar files found = %i' % len(files)

    # Remove files outside date range
    time_str = [re.search(date_str, f).group() for f in files]
    times = [datetime.strptime(string, date_fmt) for string in time_str]
    files = [
        f for f, time in zip(files, times) if time >= start and time <= stop]

    if verbose:
        print 'Number of radar files after within date range = %i' % len(files)

    # Loop over all files
    clutter = []
    for f in files:

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Read radar data
        radar = read(f, exclude_fields=exclude_fields)

        # Check radar VCP
        if vcp_sweeps is not None and radar.nsweeps != vcp_sweeps:
            continue
        if vcp_rays is not None and radar.nrays != vcp_rays:
            continue

        # Mask incoherent echoes
        ncp = radar.fields[ncp_field]['data']
        ncp = np.ma.masked_less(ncp, min_ncp, copy=False)

        # Find pixels that have a coherent signal
        clutter.append(np.logical_not(ncp.mask))

    # Compute the probability a pixel (gate) is clutter
    sample_size = len(clutter)
    clutter_map = np.sum(clutter, axis=0).astype(np.float64) / sample_size

    # Add clutter frequency map to radar object
    clutter = {
        'data': clutter_map,
        'long_name': 'Clutter frequency map',
        'standard_name': 'clutter_frequency_map',
        'valid_min': 0.0,
        'valid_max': 1.0,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field('clutter_frequency_map', clutter, replace_existing=False)

    return {
        'clutter frequency map': clutter_map,
        'last radar': radar,
        'sample size': sample_size,
        'radar files': [os.path.basename(f) for f in files],
        'min normalized coherent power': min_ncp,
        'sweeps in VCP': vcp_sweeps,
    }


def map_from_json(filename, inpdir=None, min_ncp=0.5, vcp_sweeps=None,
                  vcp_rays=None, exclude_fields=None, ncp_field=None,
                  debug=False, verbose=False):
    """
    Compute the clutter frequency (probability) map from the files listed in a
    JSON file. The listed files should define a non-precipitating time period
    where (most) echoes present must be, by definition, clutter.

    Parameters
    ----------

    Optional Parameters
    -------------------

    Returns
    -------
    """

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
    clutter = []
    for f in files:

        # Read radar data
        radar = read(f, exclude_fields=exclude_fields)

        # Check radar VCP
        if vcp_sweeps is not None and radar.nsweeps != vcp_sweeps:
            continue
        if vcp_rays is not None and radar.nrays != vcp_rays:
            continue

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Mask incoherent echoes
        ncp = radar.fields[ncp_field]['data']
        ncp = np.ma.masked_less(ncp, min_ncp, copy=False)

        # Find pixels that have a coherent signal
        clutter.append(np.logical_not(ncp.mask))

    # Compute the probability a pixel (gate) is clutter
    sample_size = len(clutter)
    clutter_map = np.sum(clutter, axis=0).astype(np.float64) / sample_size

    # Add clutter frequency map to radar object
    clutter = {
        'data': clutter_map,
        'long_name': 'Clutter frequency map',
        'standard_name': 'clutter_frequency_map',
        'valid_min': 0.0,
        'valid_max': 1.0,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field('clutter_frequency_map', clutter, replace_existing=False)

    return {
        'clutter frequency map': clutter_map,
        'last radar': radar,
        'sample size': sample_size,
        'radar files': [os.path.basename(f) for f in files],
        'min normalized coherent power': min_ncp,
        'sweeps in VCP': vcp_sweeps,
    }
