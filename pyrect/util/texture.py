"""
echo.util.texture
=================

A submodule for computing texture fields from radar fields. A texture field is
defined as the standard deviation of a radar measurement within a 1-D or 2-D
window centered around a radar gate.

"""

import time
import numpy as np

from pyart.config import get_fillvalue, get_field_name

from ._texture import compute_texture


def add_textures(
        radar, fields=None, gatefilter=None, window=(3, 3), min_sample=5,
        min_sweep=None, max_sweep=None, min_range=None, max_range=None,
        rays_wrap_around=False, fill_value=None, debug=False, verbose=False):
    """
    Add texture fields to radar.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing specified field.
    fields : str or list or tuple, optional
        Radar field(s) to compute texture field(s). If None, texture fields
        for all available radar fields will be computed and added.
    gatefilter : GateFilter, optional
        Py-ART GateFilter specifying radar gates which should be included when
        computing the texture field.
    window : list or tuple, optional
        The 2-D (ray, gate) texture window used to compute texture fields.
    min_sample : int, optional
        Minimum sample size within texture window required to define a valid
        texture. Note that a minimum of 2 radar gates are required to compute
        the texture field.
    min_sweep : int, optional
        Minimum sweep number to compute texture field.
    max_sweep : int, optional
        Maximum sweep number to compute texture field.
    min_range : float, optional
        Minimum range in meters from radar to compute texture field.
    max_range : float, optional
        Maximum range in meters from radar to compute texture field.
    fill_value : float, optional
        Value indicating missing or bad data in radar field data. If None,
        default value in Py-ART configuration file is used.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse fields to compute textures
    # If no fields are specified then the texture field of all available radar
    # fields are computed
    if fields is None:
        fields = radar.fields.keys()
    if isinstance(fields, str):
        fields = [fields]

    # Parse texture window parameters
    ray_window, gate_window = window

    if verbose:
        print 'Number of rays in window: {}'.format(ray_window)
        print 'Number of gates in window: {}'.format(gate_window)

    for field in fields:
        if verbose:
            print 'Computing texture field: {}'.format(field)

        _add_texture(
            radar, field, gatefilter=gatefilter, ray_window=ray_window,
            gate_window=gate_window, min_sample=min_sample,
            min_sweep=min_sweep, max_sweep=max_sweep, min_range=min_range,
            max_range=max_range, rays_wrap_around=rays_wrap_around,
            fill_value=fill_value, debug=debug, verbose=verbose)

    return


def _add_texture(
        radar, field, gatefilter=None, ray_window=3, gate_window=3,
        min_sample=5, min_sweep=None, max_sweep=None, min_range=None,
        max_range=None, rays_wrap_around=False, fill_value=None,
        text_field=None, debug=False, verbose=False):
    """
    Compute the texture field (standard deviation) of the input radar field
    within a 1-D or 2-D window.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing specified field.
    field : str
        Radar field to compute texture field.
    gatefilter : GateFilter, optional
        Py-ART GateFilter specifying radar gates which should be included when
        computing the texture field.
    ray_window : int, optional
        Number of rays in texture window.
    gate_window : int, optional
        Number of range gates in texture window.
    min_sample : int, optional
        Minimum sample size within texture window required to define a valid
        texture. Note that a minimum of 2 radar gates are required to compute
        the texture field.
    min_sweep : int, optional
        Minimum sweep number to compute texture field.
    max_sweep : int, optional
        Maximum sweep number to compute texture field.
    min_range : float, optional
        Minimum range in meters from radar to compute texture field.
    max_range : float, optional
        Maximum range in meters from radar to compute texture field.
    fill_value : float, optional
        Value indicating missing or bad data in radar field data. If None,
        default value in Py-ART configuration file is used.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if text_field is None:
        text_field = '{}_texture'.format(field)

    # Parse radar data
    data = radar.fields[field]['data'].copy()

    # Mask sweeps outside of specified sweep range
    for sweep, slc in enumerate(radar.iter_slice()):
        if min_sweep is not None and sweep < min_sweep:
            data[slc] = np.ma.masked
        if max_sweep is not None and sweep > max_sweep:
            data[slc] = np.ma.masked

    # Mask radar range gates outside specified gate range
    if min_range is not None:
        idx = np.abs(radar.range['data'] - min_range).argmin()
        data[:,:idx+1] = np.ma.masked
    if max_range is not None:
        idx = np.abs(radar.range['data'] - max_range).argmin()
        data[:,idx+1:] = np.ma.masked

    # Parse gate filter information
    if gatefilter is not None:
        data = np.ma.masked_where(gatefilter.gate_excluded, data)

    if debug:
        N = np.ma.count(data)
        print 'Sample size of data field: {}'.format(N)

    # Parse sweep start and end indices
    sweep_start = radar.sweep_start_ray_index['data']
    sweep_end = radar.sweep_end_ray_index['data']

    # Record starting time
    start = time.time()

    # Compute texture field
    sigma, sample_size = compute_texture(
        np.ma.filled(data, fill_value), sweep_start, sweep_end,
        ray_window=ray_window, gate_window=gate_window,
        rays_wrap_around=rays_wrap_around, fill_value=fill_value,
        debug=debug, verbose=verbose)

    # Record elapsed time to compute texture
    elapsed = time.time() - start
    if debug:
        print('Elapsed time to compute texture: {:.2f} sec'.format(elapsed))

    if min_sample is not None:
        sigma = np.ma.masked_where(sample_size < min_sample, sigma)
    sigma = np.ma.masked_invalid(sigma)
    sigma = np.ma.masked_values(sigma, fill_value, atol=1.0e-5)

    if debug:
        N = np.ma.count(sigma)
        print 'Sample size of texture field: {}'.format(N)

    sigma_dict = {
        'data': sigma,
        'units': '',
        'valid_min': 0.0,
        'number_of_rays': ray_window,
        'number_of_gates': gate_window,
        }

    radar.add_field(text_field, sigma_dict, replace_existing=True)

    return
