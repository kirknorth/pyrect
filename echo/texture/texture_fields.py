"""
clutter.texture.texture_fields
==============================

"""

import os
import json
import pickle
import numpy as np

from collections import defaultdict

from pyart.io import read
from pyart.config import get_fillvalue, get_field_name

from . import compute_texture


def _pickle_histograms(histograms, filename, outdir=None):
    """
    """

    # Parse output directory
    if outdir is not None:
        filename = os.path.join(outdir, filename)

    # Create dictionary from histogram list
    if isinstance(histograms, list):
        data = {histogram['field']: histogram for histogram in histograms}
    elif isinstance(histograms, dict):
        data = histograms
    else:
        raise ValueError('Unsupported histogram type')

    with open(filename, 'wb') as fid:
        pickle.dump(data, fid, protocol=pickle.HIGHEST_PROTOCOL)

    return


def _compute_field(radar, field, gatefilter=None, ray_window=3, gate_window=3,
                   min_sample=5, min_ncp=None, min_sweep=None, max_sweep=None,
                   min_range=None, max_range=None, rays_wrap_around=False,
                   fill_value=None, text_field=None, ncp_field=None):
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
    if text_field is None:
        text_field = '{}_texture'.format(field)
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Parse radar data
    data = radar.fields[field]['data']

    # Mask sweeps outside of specified range
    if min_sweep is not None:
        i = radar.sweep_start_ray_index['data'][min_sweep]
        data[:i+1,:] = np.ma.masked
    if max_sweep is not None:
        i = radar.sweep_end_ray_index['data'][max_sweep]
        data[i+1:,:] = np.ma.masked

    # Mask radar range gates outside specified range
    if min_range is not None:
        i = np.abs(radar.range['data'] / 1000.0 - min_range).argmin()
        data[:,:i+1] = np.ma.masked
    if max_range is not None:
        i = np.abs(radar.range['data'] / 1000.0 - max_range).argmin()
        data[:,i+1:] = np.ma.masked

    # Mask incoherent echoes
    if min_ncp is not None:
        ncp = radar.fields[ncp_field]['data']
        data = np.ma.masked_where(ncp < min_ncp, data)

    # Mask excluded gates
    if gatefilter is not None:
        data = np.ma.masked_where(gatefilter.gate_excluded, data)

    # Prepare data for ingest into Fortran wrapper
    data = np.ma.filled(data, fill_value)
    data = np.asfortranarray(data, dtype=np.float64)

    # Parse sweep start/end indices
    # Offset indices in order to be compatible with Fortran and avoid a
    # segmentation fault
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
    texture = np.ma.masked_invalid(texture, copy=False)
    texture.set_fill_value(fill_value)

    # Create texture field dictionary and add it to the radar object
    texture = {
        'data': texture.astype(np.float32),
        'long_name': '{} texture'.format(radar.fields[field]['long_name']),
        'standard_name': text_field,
        'valid_min': 0.0,
        '_FillValue': texture.fill_value,
        'units': radar.fields[field]['units'],
        'comment_1': ('Texture field is defined as the standard deviation '
                      'within a prescribed 2D window'),
        'comment_2': '{} x {} window'.format(gate_window, ray_window),
    }
    radar.add_field(text_field, texture, replace_existing=True)

    return


def add_textures(radar, fields=None, gatefilter=None, texture_window=(3, 3),
                 texture_sample=5, min_ncp=None, min_sweep=None,
                 max_sweep=None, min_range=None, max_range=None,
                 rays_wrap_around=False, fill_value=None, ncp_field=None):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    if fields is None:
        fields = radar.fields.keys()

    # Parse texture window parameters
    ray_window, gate_window = texture_window

    for field in fields:
        _compute_field(
            radar, field, gatefilter=gatefilter, ray_window=ray_window,
            gate_window=gate_window, min_sample=texture_sample,
            min_ncp=min_ncp, min_sweep=min_sweep, max_sweep=max_sweep,
            min_range=min_range, max_range=max_range,
            rays_wrap_around=rays_wrap_around, fill_value=fill_value,
            text_field=None, ncp_field=ncp_field)

    return


def histogram_from_json(
        filename, field, inpdir=None, texture_window=(3, 3), min_sample=5,
        num_bins=10, limits=None, min_ncp=0.5, vcp_sweeps=None,
        vcp_rays=None, min_sweep=None, max_sweep=None, min_range=None,
        max_range=None, rays_wrap_around=False, exclude_fields=None,
        fill_value=None, ncp_field=None, verbose=False):
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

    # Parse texture window parameters
    ray_window, gate_window = texture_window

    # Loop over all files
    counts = np.zeros(num_bins, dtype=np.float64)
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

        # Compute texture fields
        _compute_field(
            radar, field, ray_window=ray_window, gate_window=gate_window,
            min_sample=min_sample, min_ncp=min_ncp, min_sweep=min_sweep,
            max_sweep=max_sweep, min_range=min_range, max_range=max_range,
            rays_wrap_around=rays_wrap_around, fill_value=fill_value,
            ncp_field=ncp_field)

        # Parse data and compute histogram
        data = radar.fields['{}_texture'.format(field)]['data']
        hist, bin_edges = np.histogram(
            data.compressed(), bins=num_bins, range=limits, normed=False,
            weights=None, density=False)
        counts += hist

    # Compute bin centers
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0

    # Compute normalized histogram and probability density
    counts_norm = counts / counts.max()
    pdf = counts_norm / np.sum(counts_norm * np.diff(bin_edges))

    return {
        'field': '{}_texture'.format(field),
        'histogram counts': counts,
        'normalized histogram': counts_norm,
        'probability density': pdf,
        'number of bins': num_bins,
        'limits': limits,
        'bin edges': bin_edges,
        'bin centers': bin_centers,
        'radar files': [os.path.basename(f) for f in files],
        'min sweep': min_sweep,
        'max sweep': max_sweep,
        'min range': min_range,
        'max range': max_range,
        'min normalized coherent power': min_ncp,
        'sweeps in VCP': vcp_sweeps,
        'rays in VCP': vcp_rays,
        'ray window size': ray_window,
        'gate window size': gate_window,
        }

def histograms_from_radar(
        radar, hist_dict, gatefilter=None, texture_window=(3, 3),
        texture_sample=5, min_ncp=None, min_sweep=None, max_sweep=None,
        min_range=None, max_range=None, rays_wrap_around=False,
        fill_value=None, ncp_field=None, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # TODO: check input histogram dictionary for proper keys

    # Parse texture window parameters
    ray_window, gate_window = texture_window

    # Loop over all fields and compute histogram counts
    for field in hist_dict:

        # Compute texture fields
        _compute_field(
            radar, field, gatefilter=gatefilter, ray_window=ray_window,
            gate_window=gate_window, min_sample=texture_sample,
            min_ncp=min_ncp, min_sweep=min_sweep, max_sweep=max_sweep,
            min_range=min_range, max_range=max_range,
            rays_wrap_around=rays_wrap_around, fill_value=fill_value,
            ncp_field=ncp_field)

        # Parse histogram parameters
        bins = hist_dict[field]['number of bins']
        limits = hist_dict[field]['limits']

        # Parse data and compute histogram
        data = radar.fields['{}_texture'.format(field)]['data']
        counts, bin_edges = np.histogram(
            data.compressed(), bins=bins, range=limits, normed=False,
            weights=None, density=False)
        hist_dict[field]['histogram counts'] += counts

        # Compute bin centers and add to dictionary
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
        hist_dict[field]['bin edges'] = bin_edges
        hist_dict[field]['bin centers'] = bin_centers

    return
