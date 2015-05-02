"""
echo.correct.basic_fixes
========================

"""

import numpy as np

from pyart.config import get_fillvalue, get_field_name

from . import sweeps


def remove_salt(radar, fields=None, salt_window=(3, 3), salt_sample=5,
                rays_wrap_around=False, mask_data=True, fill_value=None,
                debug=False, verbose=False):
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

        if debug:
            print 'Removing salt and pepper noise from: {}'.format(field)

        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype_orig = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        if debug:
            n = np.sum(data != fill_value)
            print 'Sample size before salt removal: {}'.format(n)

        # Fortran wrapper
        sweeps.remove_salt(
            data, sweep_start, sweep_end, ray_window=ray_window,
            gate_window=gate_window, min_sample=salt_sample,
            rays_wrap=rays_wrap_around, fill_value=fill_value)

        if debug:
            n = np.sum(data != fill_value)
            print 'Sample size after salt removal: {}'.format(n)

        # Mask invalid data
        if mask_data:
            data = np.ma.masked_equal(data, fill_value, copy=False)

        radar.fields[field]['data'] = data.astype(dtype_orig)

    return


def interpolate_missing(
        radar, fields=None, ray_window=3, gate_window=3, min_sample=8,
        kind='mean', rays_wrap_around=False, fill_value=None, debug=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names to interpolate
    if fields is None:
        fields = radar.fields.keys()

    for field in fields:

        if debug:
            print 'Filling missing gates for field: {}'.format(field)

        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype_orig = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        if debug:
            n = np.sum(data == fill_value)
            print 'Number of missing gates before fill : {}'.format(n)

        # Parse sweep parameters
        # Offset index arrays in order to be compatible with Fortran and avoid
        # a segmentation fault
        sweep_start = radar.sweep_start_ray_index['data'] + 1
        sweep_end = radar.sweep_end_ray_index['data'] + 1

        # Call Fortran routine
        if kind.upper() == 'MEAN':
            sweeps.mean_fill(
                data, sweep_start, sweep_end, ray_window=ray_window,
                gate_window=gate_window, min_sample=min_sample,
                rays_wrap=rays_wrap_around, fill_value=fill_value)
        else:
            raise ValueError('Unsupported interpolation method')

        if debug:
            n = np.sum(data == fill_value)
            print 'Number of missing gates after fill : {}'.format(n)

        # Mask invalid data
        data = np.ma.masked_equal(data, fill_value, copy=False)

        # Add interpolated data to radar object
        radar.fields[field]['data'] = data.astype(dtype_orig)

    return


def _binary_dilation(radar, field, structure=None, iterations=1):
    """
    """

    for sweep in radar.iter_slice():

        # Parse radar data
        mask = radar.fields[field]['data'][sweep]
        mask = np.ma.filled(mask, False)

        # Call SciPy's binary fill algorithm
        mask = ndimage.binary_dilation(
            mask, structure=structure, iterations=iterations, mask=None,
            output=None, border_value=0, origin=0, brute_force=False)

        # Update radar sweep
        radar.fields[field]['data'][sweep] = mask

    return


def _binary_fill(radar, field, structure=None):
    """
    """

    # Parse radar data
    mask = radar.fields[field]['data']
    mask = np.ma.filled(mask, False)

    # Call SciPy's binary fill algorithm
    mask = ndimage.binary_fill_holes(
        mask, structure=structure, output=None, origin=0)

    # Update radar field
    radar.fields[field]['data'] = mask

    return
