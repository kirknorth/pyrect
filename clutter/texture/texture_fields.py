"""
clutter.texture.texture_fields
==============================

"""

import numpy as np

from pyart.config import get_fillvalue

from . import compute_texture


def compute_field(radar, field, ray_window=3, gate_window=3, min_sample=None,
                  fill_value=None):
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

    # Parse data field and sweep parameters
    # Need to add 1 to sweep index arrays in order to be consistent with
    # Fortran array indexing
    data = np.ma.filled(radar.fields[field]['data'], fill_value)
    sweep_start = radar.sweep_start_ray_index['data'] + 1
    sweep_end = radar.sweep_end_ray_index['data'] + 1

    sample_size, texture = compute_texture.compute(
        data, sweep_start, sweep_end, ray_window=ray_window,
        gate_window=gate_window, fill_value=fill_value)

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
                    'of a field within a prescribed window')
        'comment_2': 'Window size was {} x {}'.format(gate_window, ray_window)
    }

    radar.add_field('{}_texture'.format(field), texture, replace_existing=True)

    return
