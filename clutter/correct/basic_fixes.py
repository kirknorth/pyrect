"""
clutter.correct.basic_fixes
===========================

"""

import numpy as np

from pyart.config import get_fillvalue

from . import fill_holes

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

        # Parse sweep parameters
        # Need to add 1 to sweep index arrays in order to be consistent with
        # Fortran array indexing and avoid a segmentation fault
        sweep_start = radar.sweep_start_ray_index['data'] + 1
        sweep_end = radar.sweep_end_ray_index['data'] + 1

        # Parse radar data
        data = radar.fields[field]['data']
        data = np.ma.filled(data, fill_value).astype(np.float64)

        # Call Fortran routine
        if kind == 'mean':
            fill_holes.mean_fill(
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
