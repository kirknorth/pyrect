"""
clutter.correct.basic_fixes
===========================

"""

import numpy as np

from pyart.config import get_fillvalue, get_field_name

from . import sweeps


def remove_salt(radar, fields=None, salt_window=(3, 3), salt_sample=5,
                rays_wrap_around=False, mask_data=True, fill_value=None,
                verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if fields is None:
        fields = radar.fields.keys()

    # Parse sweep start/end indices
    # Offset indices in order to be compatible with Fortran and avoid a
    # segmentation fault
    sweep_start = radar.sweep_start_ray_index['data'] + 1
    sweep_end = radar.sweep_end_ray_index['data'] + 1

    # Parse window size
    ray_window, gate_window = salt_window

    # Remove salt and pepper noise for each field
    for field in fields:
        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        if verbose:
            sample = np.sum(data != fill_value)
            print 'Sample size before salt removal = %i' % sample

        # Fortran wrapper
        sweeps.remove_salt(
            data, sweep_start, sweep_end, ray_window=ray_window,
            gate_window=gate_window, min_sample=salt_sample,
            rays_wrap=rays_wrap_around, fill_value=fill_value)

        if verbose:
            sample = np.sum(data != fill_value)
            print 'Sample size after salt removal = %i' % sample

        # Mask invalid data
        if mask_data:
            data = np.ma.masked_equal(data, fill_value, copy=False)

        radar.fields[field]['data'] = data.astype(dtype)

    return


def interpolate_missing(radar, fields=None, ray_window=3, gate_window=3,
                        min_sample=8, kind='mean', fill_value=None):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names to interpolate
    if fields is None:
        fields = radar.fields.keys()

    for field in fields:

        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        # Parse sweep parameters
        # Offset index arrays in order to be compatible with Fortran and avoid
        # a segmentation fault
        sweep_start = radar.sweep_start_ray_index['data'] + 1
        sweep_end = radar.sweep_end_ray_index['data'] + 1

        # Call Fortran routine
        if kind == 'mean':
            sweeps.mean_fill(
                data, sweep_start, sweep_end, ray_window=ray_window,
                gate_window=gate_window, min_sample=min_sample,
                fill_value=fill_value)
        else:
            raise ValueError('Unsupported interpolation method')

        # Mask invalid data
        data = np.ma.masked_equal(data, fill_value, copy=False)

        # Add interpolated data to radar object
        radar.fields[field]['data'] = data

    return
