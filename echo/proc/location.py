"""
echo.location.geo
=================

"""

import os
import pickle
import numpy as np

from pyart.config import get_fillvalue, get_field_name

from . import common


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


def height_histogram_from_radar(
        radar, hist_dict, gatefilter=None, min_ncp=None, min_sweep=None,
        max_sweep=None, min_range=None, max_range=None, fill_value=None,
        ncp_field=None, verbose=False, debug=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # TODO: check input histogram dictionary for proper keys

    # Compute locations of radar gates and convert to kilometers
    x, y, heights = common.standard_refraction(radar, use_km=True)

    if debug:
        print '(Min, Max) radar gate height = ({:.3f}, {:.3f}) km'.format(
            heights.min(), heights.max())

    # Mask sweeps outside specified range
    if min_sweep is not None:
        i = radar.sweep_start_ray_index['data'][min_sweep]
        heights[:i+1,:] = np.ma.masked
    if max_sweep is not None:
        i = radar.sweep_end_ray_index['data'][max_sweep]
        heights[i+1:,:] = np.ma.masked

    # Mask radar range gates outside specified range
    if min_range is not None:
        i = np.abs(radar.range['data'] / 1000.0 - min_range).argmin()
        heights[:,:i+1] = np.ma.masked
    if max_range is not None:
        i = np.abs(radar.range['data'] / 1000.0 - max_range).argmin()
        heights[:,i+1:] = np.ma.masked

    # Mask incoherent echoes
    if min_ncp is not None:
        ncp = radar.fields[ncp_field]['data']
        heights = np.ma.masked_where(ncp < min_ncp, heights)

    # Mask excluded gates
    if gatefilter is not None:
        heights = np.ma.masked_where(gatefilter.gate_excluded, heights)

    # Parse histogram parameters
    bins = hist_dict['number of bins']
    limits = hist_dict['limits']

    # Compute histogram counts
    counts, bin_edges = np.histogram(
        heights.compressed(), bins=bins, range=limits, normed=False,
        weights=None, density=False)
    hist_dict['histogram counts'] += counts

    # Parse bin edges and bin centers
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
    hist_dict['bin edges'] = bin_edges
    hist_dict['bin centers'] = bin_centers

    return
