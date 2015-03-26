"""
clutter.moment.moment_fields
============================

"""

import os
import json
import pickle
import numpy as np

from pyart.io import read
from pyart.config import get_fillvalue, get_field_name

from ..correct import noise


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


def histogram_from_json(
        filename, field, inpdir=None, bins=10, limits=None, min_ncp=0.5,
        vcp_sweeps=None, vcp_rays=None, min_sweep=None, max_sweep=None,
        exclude_fields=None, fill_value=None, ncp_field=None, verbose=False):
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

        # Check radar VCP
        if vcp_sweeps is not None and radar.nsweeps != vcp_sweeps:
            continue
        if vcp_rays is not None and radar.nrays != vcp_rays:
            continue

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Parse radar fields
        data = radar.fields[field]['data']

        # Mask sweeps outside specified range
        if min_sweep is not None:
            i = radar.sweep_start_ray_index['data'][min_sweep]
            data[:i+1,:] = np.ma.masked
        if max_sweep is not None:
            i = radar.sweep_end_ray_index['data'][max_sweep]
            data[i+1:,:] = np.ma.masked

        # Mask incoherent echoes
        if min_ncp is not None:
            ncp = radar.fields[ncp_field]['data']
            data = np.ma.masked_where(ncp < min_ncp, data)

        # Bin data and compute frequencies
        hist, bin_edges = np.histogram(
            data.compressed(), bins=bins, range=limits, normed=False,
            weights=None, density=False)
        histogram += hist

    # Compute bin centers
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0

    # Compute normalized histogram and probability density
    histogram_norm = histogram / histogram.max()
    pdf = histogram_norm / np.sum(histogram_norm * np.diff(bin_edges))

    return {
        'field': field,
        'histogram counts': histogram,
        'normalized histogram': histogram_norm,
        'probability density': pdf,
        'number of bins': bins,
        'limits': limits,
        'bin edges': bin_edges,
        'bin centers': bin_centers,
        'radar files': [os.path.basename(f) for f in files],
        'min sweep': min_sweep,
        'max sweep': max_sweep,
        'min normalized coherent power': min_ncp,
        'sweeps in VCP': vcp_sweeps,
        'rays in VCP': vcp_rays,
        }


def histograms_from_radar(
        radar, hist_dict, gatefilter=None, min_ncp=None, min_sweep=None,
        max_sweep=None, min_range=None, max_range=None, fill_value=None,
        ncp_field=None, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # TODO: check input histogram dictionary for proper keys

    # Loop over all fields and compute histogram counts
    for field in hist_dict:

        # Parse radar fields
        data = radar.fields[field]['data']

        # Mask sweeps outside specified range
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
            data = np.ma.masked_where(
                gatefilter.gate_excluded, data)

        # Parse histogram parameters
        bins = hist_dict[field]['number of bins']
        limits = hist_dict[field]['limits']

        # Bin data and compute frequencies
        counts, bin_edges = np.histogram(
            data.compressed(), bins=bins, range=limits, normed=False,
            weights=None, density=False)
        hist_dict[field]['histogram counts'] += counts

        # Parse bin edges and bin centers
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
        hist_dict[field]['bin edges'] = bin_edges
        hist_dict[field]['bin centers'] = bin_centers

    return
