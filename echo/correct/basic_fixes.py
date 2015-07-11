"""
echo.correct.basic_fixes
========================

"""

import numpy as np

from scipy import ndimage

from pyart.config import get_fillvalue, get_field_name

from . import sweeps


def remove_salt(radar, fields=None, salt_window=(3, 3), salt_sample=5,
                rays_wrap_around=False, mask_data=True, fill_value=None,
                debug=False, verbose=False):
    """
    Remove basic salt and pepper noise from radar fields. Noise removal is done
    in-place for each radar field, i.e., no new field is created, rather
    the original field is changed.

    Parameters
    ----------
    radar : Radar
        Radar object containing the specified fields for noise removal.
    fields : str, list or tuple, optional
        The field(s) which will have basic salt and pepper noised removed.
    salt_window : tuple or list, optional
        The 2-D (ray, gate) window filter used to determine whether a radar
        gate is isolated (noise) or part of a larger feature.
    salt_sample : int, optional
        The minimum sample size within 'salt_window' for a radar gate to be
        considered part of a larger feature. If the sample size within
        salt_window is below this value, then the radar gate is considered to
        be isolated and therefore salt and pepper noise.
    rays_wrap_around : bool, optional

    mask_data : bool, optional
        Whether the radar field(s) should be masked after salt and pepper noise
        is removed. This should be set to False when the field in question is
        a binary (mask) field, e.g., radar significant detection mask.
    fill_value : float, optional
        The fill value for radar fields. If not specified, the default fill
        value from the Py-ART configuration is used.
    debug, verbose : bool, optional
        True to print debugging and progress information, respectively, False
        to suppress.

    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if fields is None:
        fields = radar.fields.keys()

    # Check if input fields is single field
    if isinstance(fields, str):
        fields = [fields]

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
            print 'Removing salt and pepper noise: {}'.format(field)

        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype_orig = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        if debug:
            N = np.count_nonzero(~np.isclose(data, fill_value, atol=1.0e-5))
            print 'Sample size before salt removal: {}'.format(N)

        # Fortran wrapper
        sweeps.remove_salt(
            data, sweep_start, sweep_end, ray_window=ray_window,
            gate_window=gate_window, min_sample=salt_sample,
            rays_wrap=rays_wrap_around, fill_value=fill_value)

        if debug:
            N = np.count_nonzero(~np.isclose(data, fill_value, atol=1.0e-5))
            print 'Sample size after salt removal: {}'.format(N)

        # Mask invalid data
        if mask_data:
            data = np.ma.masked_equal(data, fill_value, copy=False)
            data.set_fill_value(fill_value)

        radar.fields[field]['data'] = data.astype(dtype_orig)

    return


def interpolate_missing(
        radar, fields=None, interp_window=(3, 3), interp_sample=8,
        kind='mean', rays_wrap_around=False, fill_value=None, debug=False,
        verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names to interpolate
    if fields is None:
        fields = radar.fields.keys()

    # Parse interpolation parameters
    ray_window, gate_window = interp_window

    # Loop over all fields and interpolate missing gates
    for field in fields:

        if verbose:
            print 'Filling missing gates: {}'.format(field)

        # Parse radar data and its original data type
        data = radar.fields[field]['data']
        dtype_orig = data.dtype

        # Prepare data for ingest into Fortran wrapper
        data = np.ma.filled(data, fill_value)
        data = np.asfortranarray(data, dtype=np.float64)

        if debug:
            N = np.count_nonzero(~np.isclose(data, fill_value, atol=1.0e-5))
            print 'Sample size before fill: {}'.format(N)

        # Parse sweep parameters
        # Offset index arrays in order to be compatible with Fortran and avoid
        # a segmentation fault
        sweep_start = radar.sweep_start_ray_index['data'] + 1
        sweep_end = radar.sweep_end_ray_index['data'] + 1

        # Call Fortran routine
        if kind.upper() == 'MEAN':
            sweeps.mean_fill(
                data, sweep_start, sweep_end, ray_window=ray_window,
                gate_window=gate_window, min_sample=interp_sample,
                rays_wrap=rays_wrap_around, fill_value=fill_value)
        else:
            raise ValueError('Unsupported interpolation method')

        if debug:
            N = np.count_nonzero(~np.isclose(data, fill_value, atol=1.0e-5))
            print 'Sample size after fill: {}'.format(N)

        # Mask invalid data
        data = np.ma.masked_equal(data, fill_value, copy=False)
        data.set_fill_value(fill_value)

        # Add interpolated data to radar object
        radar.fields[field]['data'] = data.astype(dtype_orig)

    return


def _binary_dilation(radar, field, structure=None, iterations=1, debug=False,
                     verbose=False):
    """
    Use binary dilation to expand around the edges of a binary (mask) radar
    field. Non-zero (True) elements of the radar field form the subset to be
    dilated. Unexpected results may occur if the specified field is not a
    binary field.

    Parameters
    ----------
    radar : Radar
        Radar object containing the specified field to be dilated.
    field : str
        Radar field name to be dilated.
    structure : array_like, optional
        Binary structuring element used to dilate radar field. The default
        structuring element has a squared connectivity equal to one.
    iterations : int, optional
        The number of iterations to repeat dilation. If iterations is less than
        1, the dilation is repeated until the result does not change anymore.
    debug, verbose : bool, optional
        True to print debugging and progress information, respectively, False
        to suppress.

    """

    if verbose:
        print 'Performing binary dilation: {}'.format(field)

    # Parse binary structuring element
    if structure is None:
        structure = ndimage.generate_binary_structure(2, 1)

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size before dilation: {}'.format(N)

    for sweep in radar.iter_slice():

        # Parse radar sweep data
        # Non-zero elements of the array form the subset to be dilated
        is_valid = np.ma.filled(radar.fields[field]['data'][sweep], 0)

        # Call SciPy's binary dilation algorithm
        is_valid = ndimage.binary_dilation(
            is_valid, structure=structure, iterations=iterations, mask=None,
            output=None, border_value=0, origin=0, brute_force=False)

        radar.fields[field]['data'][sweep] = is_valid

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size after dilation: {}'.format(N)

    return


def _binary_significant_features(
        radar, binary_field, size_bins=75, size_limits=(0, 300),
        structure=None, debug=False, verbose=False):
    """
    Objectively determine the minimum echo feature size (in radar gates) and
    remove features smaller than this from the specified radar field. This
    function can be used to objectively remove salt and pepper noise from
    binary (mask) radar fields. Unexpected results may occur if the specified
    radar field is not a binary field.

    Parameters
    ----------
    radar : Radar
        Radar object containing the specified binary field.
    binary_field : str
        The binary radar field that will have insignificant echo features
        removed.
    size_bins : int, optional
        Number of bins used to bin echo feature sizes and thus define its
        distribution.
    size_limits : list or tuple, optional
        Limits of the echo feature size distribution. The upper limit needs to
        be large enough to include the minimum feature size. The size bin width
        is defined by both the size_bins and size_limits parameters.
    structure : array_like, optional
        Binary structuring element used to define connected features. The
        default structuring element has a squared connectivity equal to one.
    debug, verbose : bool, optional
        True to print debugging and progress information, respectively, False
        to suppress.

    """

    # Parse binary structuring element
    if structure is None:
        structure = ndimage.generate_binary_structure(2, 1)

    # Parse feature size arrays
    size_data = np.zeros_like(
        radar.fields[binary_field]['data'], dtype=np.int32)
    feature_sizes = []

    for i, sweep in enumerate(radar.iter_slice()):

        # Parse radar sweep data
        # Non-zero elements of the array form the subset to be dilated
        data = radar.fields[binary_field]['data'][sweep]
        is_valid_gate = np.ma.filled(data, 0)

        # Label the connected features in the sweep data and create index
        # array which defines each unique labeled feature
        labels, nlabels = ndimage.label(
            is_valid_gate, structure=structure, output=None)
        index = np.arange(1, nlabels + 1, 1)

        if debug:
            print 'Unique features in sweep {}: {}'.format(i, nlabels)

        # Compute the size in radar gates of each labeled feature
        sweep_sizes = ndimage.labeled_comprehension(
            is_valid_gate, labels, index, np.count_nonzero, np.int32, 0)
        feature_sizes.append(sweep_sizes)

        # Set each labeled feature to its total size in radar gates
        for label, size in zip(index, sweep_sizes):
            size_data[sweep][labels == label] = size

    feature_sizes = np.hstack(feature_sizes)

    # Bin and count occurrences of labeled feature sizes
    # Compute bin centers and bin width
    counts, bin_edges = np.histogram(
        feature_sizes, bins=size_bins, range=size_limits, normed=False,
        weights=None, density=False)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
    bin_width = np.diff(bin_edges).mean()

    if debug:
        print 'Feature size bin width: {} gate(s)'.format(bin_width)

    # Compute the peak of the labeled feature size distribution
    # We expect the peak of this distribution to be close to 1 radar gate
    peak_size = bin_centers[counts.argmax()] - bin_width / 2.0

    if debug:
        print 'Feature size at peak: {} gate(s)'.format(peak_size)

    # Determine the first instance when the count (sample size) of a feature
    # size bin reaches 0 in the right side of the feature size distribution
    # This will define the minimum feature size
    is_zero_size = np.logical_and(
        bin_centers > peak_size, np.isclose(counts, 0, atol=1.0e-5))
    min_size = bin_centers[is_zero_size].min() - bin_width / 2.0

    if debug:
        _range = [0.0, min_size]
        print 'Insignificant feature size range: {} gates'.format(_range)

    # Remove insignificant features from the binary radar field
    radar.fields[binary_field]['data'][size_data < min_size] = 0

    return


def _binary_fill(radar, field, structure=None, debug=False, verbose=False):
    """
    """

    # Parse binary structuring element
    if structure is None:
        structure = ndimage.generate_binary_structure(2, 1)

    for sweep in radar.iter_slice():

        is_valid = np.ma.filled(radar.fields[field]['data'][sweep], 0)

        # Call SciPy's binary fill algorithm
        is_valid = ndimage.binary_fill_holes(
            is_valid, structure=structure, output=None, origin=0)

        radar.fields[field]['data'][sweep] = is_valid

    return
