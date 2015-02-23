"""
clutter.texture.texture_fields
==============================

"""

import os
import json
import numpy as np

from collections import defaultdict

from pyart.config import get_fillvalue, get_field_name
from pyart.io import read

from . import compute_texture


def _pickle_texture_histograms(histograms, filename, outdir=None):
    """
    """

    # Parse output directory
    if outdir is not None:
        filename = os.path.join(outdir, filename)

    # Create dictionary containing all texture fields organized by field name
    data = {}
    for histogram in histograms:
        data[histogram['field']] = histogram

    with open(filename, 'wb') as fid:
        pickle.dump(data, fid, protocol=pickle.HIGHEST_PROTOCOL)

    return


def _compute_field(radar, field, ray_window=3, gate_window=3, min_sample=None,
                   min_ncp=0.3, fill_value=None, ncp_field=None):
    """
    Compute the texture (standard deviation) within the 2-D window for the
    specified field.

    Parameters
    ----------

    Optional Parameters
    ----------------

    Returns
    -------
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Parse fields
    ncp = radar.fields[ncp_field]['data']
    data = radar.fields[field]['data']

    # Mask incoherent echoes
    data = np.ma.masked_where(ncp < min_ncp, data, copy=False)
    data = np.ma.filled(data, fill_value)

    # Parse sweep parameters
    # Need to add 1 to sweep index arrays in order to be consistent with
    # Fortran array indexing
    sweep_start = radar.sweep_start_ray_index['data'] + 1
    sweep_end = radar.sweep_end_ray_index['data'] + 1

    # Compute texture field
    sample_size, texture = compute_texture.compute(
        data, sweep_start, sweep_end, ray_window=ray_window,
        gate_window=gate_window, fill_value=fill_value)

    # Mask pixels (gates) where the sample size used to compute the texture
    # field was too small
    if min_sample is not None:
        texture = np.ma.masked_where(
            sample_size < min_sample, texture, copy=False)

    # Mask invalid values
    texture = np.ma.masked_equal(texture, fill_value, copy=False)

    texture = {
        'data': texture,
        'standard_name': '{}_texture'.format(
            radar.fields[field]['standard_name']),
        'long_name': '{} texture'.format(radar.fields[field]['long_name']),
        '_FillValue': texture.fill_value,
        'units': radar.fields[field]['units'],
        'comment_1': ('Texture field is defined as the standard deviation '
                    'of a field within a prescribed window'),
        'comment_2': 'Window size was {} x {}'.format(gate_window, ray_window),
    }

    radar.add_field('{}_texture'.format(field), texture, replace_existing=True)

    return


def add_textures(radar, fields, ray_window=3, gate_window=3, min_sample=None,
                 min_ncp=0.3, min_sweep=None, max_sweep=None, fill_value=None,
                 ncp_field=None):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_field_name()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    if min_sweep is not None:
        i = radar.sweep_start_ray_index['data'][min_sweep]
        radar.fields[field]['data'][0:i,:] = np.ma.masked

    if max_sweep is not None:
        i = radar.sweep_end_ray_index['data'][max_sweep]
        radar.fields[field]['data'][i+1:-1,:] = np.ma.masked

    for field in fields:
        _compute_field(
            radar, field, ray_window=ray_window, gate_window=gate_window,
            min_sample=min_sample, min_ncp=min_ncp, fill_value=fill_value,
            ncp_field=ncp_field)

    return


def histogram_from_json(
        filename, field, inpdir=None, ray_window=3, gate_window=3,
        min_sample=None, bins=10, limits=None, min_ncp=0.3, vcp_sweeps=22,
        min_sweep=None, max_sweep=None, exclude_fields=None, fill_value=None,
        ncp_field=None, verbose=False):
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

        if min_sweep is not None:
            i = radar.sweep_start_ray_index['data'][min_sweep]
            radar.fields[field]['data'][0:i,:] = np.ma.masked

        if max_sweep is not None:
            i = radar.sweep_end_ray_index['data'][max_sweep]
            radar.fields[field]['data'][i+1:-1,:] = np.ma.masked

        # Compute texture fields
        _compute_field(
            radar, field, ray_window=ray_window, gate_window=gate_window,
            min_sample=min_sample, min_ncp=min_ncp, fill_value=fill_value,
            ncp_field=ncp_field)

        # Parse data and compute histogram
        data = radar.fields['{}_texture'.format(field)]['data']
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
        }
