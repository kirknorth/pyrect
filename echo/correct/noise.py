"""
echo.correct.noise
==================

"""

import numpy as np

from scipy import ndimage
from scipy.signal import argrelextrema

from pyart.correct import GateFilter
from pyart.config import get_fillvalue, get_field_name, get_metadata

from . import sweeps, basic_fixes
from ..texture import texture_fields


def significant_detection(
        radar, gatefilter=None, remove_small_features=True, size_bins=75,
        size_limits=(0, 300), fill_holes=False, dilate=False, structure=None,
        iterations=1, rays_wrap_around=False, min_ncp=None, ncp_field=None,
        detect_field=None, debug=False, verbose=False):
    """
    Determine the significant detection of a radar. Note that significant
    detection can still include other non-meteorological echoes that the user
    may still have to remove further down the processing chain.

    Parameters
    ----------
    radar : Radar
        Radar object used to determine the appropriate GateFilter.
    gatefilter : GateFilter, optional
        If None, all radar gates will initially be assumed valid.
    remove_small_features : bool, optional
        True to remove insignificant echo features (e.g., salt and pepper
        noise) from significant detection mask.
    size_bins : int, optional
        Number of bins used to bin echo feature sizes and thus define its
        distribution.
    size_limits : list or tuple, optional
        Limits of the echo feature size distribution. The upper limit needs to
        be large enough to include the minimum feature size.
    fill_holes : bool, optional
        Fill any holes in the significant detection mask. For most radar
        volumes this should not be used since the default structuring element
        will automatically fill any sized hole.
    dilate : bool, optional
        Use binary dilation to fill in edges of the significant detection mask.
    structure : array_like, optional
        The binary structuring element used for all morphology routines. See
        SciPy's ndimage documentation for more information.
    iterations : int, optional
        The number of iterations to repeat binary dilation. If iterations is
        less than 1, binary dilation is repeated until the result does not
        change anymore.
    rays_wrap_around : bool, optional
        Whether the rays at the beginning and end of a sweep are connected
        (e.g., PPI VCP).
    min_ncp : float, optional
        Minimum normalized coherent power (signal quality) value used to
        indicate a significant echo.
    ncp_field : str, optional
        Minimum normalized coherent power (signal quality) field name. The
        default uses the Py-ART configuation file.
    detect_field : str, optional
        Radar significant detection mask field name.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress information, False to suppress.

    Returns
    -------
    gatefilter : GateFilter
        Py-ART GateFilter object indicating which radar gates are valid and
        invalid.
    """

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if detect_field is None:
        detect_field = 'significant_detection_mask'

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Exclude gates with poor signal quality
    if min_ncp is not None and ncp_field in radar.fields:
        gatefilter.include_above(ncp_field, min_ncp, op='and', inclusive=True)

    detect_dict = {
        'data': gatefilter.gate_included.astype(np.int8),
        'long_name': 'Radar significant detection mask',
        'standard_name': 'significant_detection_mask',
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = no significant detection, 1 = significant detection',
    }
    radar.add_field(detect_field, detect_dict, replace_existing=True)

    # Remove insignificant features from significant detection mask
    if remove_small_features:
        basic_fixes._binary_significant_features(
            radar, detect_field, size_bins=size_bins, size_limits=size_limits,
            structure=structure, debug=debug, verbose=verbose)

    # Fill holes in significant detection mask
    if fill_holes:
        basic_fixes._binary_fill(radar, detect_field, structure=structure)

    # Dilate significant detection mask
    if dilate:
        basic_fixes._binary_dilation(
            radar, detect_field, structure=structure, iterations=iterations,
            debug=debug, verbose=verbose)

    # Update gate filter
    gatefilter.include_equal(detect_field, 1, op='new')

    return gatefilter


def hildebrand_noise(
        radar, gatefilter=None, scale=1.0, remove_small_features=False,
        size_bins=75, size_limits=(0, 300), rays_wrap_around=False,
        fill_value=None, power_field=None, noise_field=None, mask_field=None,
        verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if power_field is None:
        power_field = get_field_name('signal_to_noise_ratio')
    if noise_field is None:
        noise_field = 'radar_noise_floor'
    if mask_field is None:
        mask_field = 'radar_noise_floor_mask'

    # Parse radar power data
    power = radar.fields[power_field]['data']

    # Prepare data for ingest into Fortran wrapper
    power = np.ma.filled(power, fill_value)
    power = np.asfortranarray(power, dtype=np.float64)

    # Convert power units to linear and sort in descending order
    # TODO: check if units are already linear
    power = np.where(power != fill_value, 10.0**(power / 10.0), power)
    power = np.sort(power, axis=0, kind='mergesort')[::-1]

    # Fortran wrapper to Hildebrand and Sekhon (1974) algorithm
    P, Q, R2, N = sweeps.hildebrand(power, fill_value=fill_value)

    # Mask invalid data
    P = np.ma.masked_equal(P, fill_value)
    Q = np.ma.masked_equal(Q, fill_value)
    R2 = np.ma.masked_equal(R2, fill_value)
    N = np.ma.masked_equal(N, fill_value)

    # Estimate noise floor in decibels and tile to proper dimensions
    noise = 10.0 * np.ma.log10(P + scale * np.ma.sqrt(Q))
    noise = np.tile(noise, (radar.nrays, 1))
    noise.set_fill_value(fill_value)

    # Add Hildebrand noise floor field to radar
    noise_dict = {
        'data': noise.astype(np.float32),
        'long_name': 'Radar noise floor estimate',
        'standard_name': noise_field,
        'units': 'dB',
        '_FillValue': noise.fill_value,
        'comment': ('Noise floor is estimated using Hildebrand and '
                    'Sekhon (1974) algorithm'),
    }
    radar.add_field(noise_field, noise_dict, replace_existing=True)

    # Compute noise floor mask
    power = radar.fields[power_field]['data']
    noise = radar.fields[noise_field]['data']
    is_noise = np.ma.filled(power >= noise, False)

    # Create radar noise floor mask dictionary
    mask_dict = {
        'data': is_noise.astype(np.int8),
        'long_name': 'Noise floor mask',
        'standard_name': mask_field,
        'valid_min': 0,
        'valid_max': 1,
        'units': 'unitless',
        '_FillValue': None,
        'comment': '0 = below noise floor, 1 = at or above noise floor',
    }
    radar.add_field(mask_field, mask_dict, replace_existing=True)

    # Remove insignificant features from noise floor mask
    if remove_small_features:
        basic_fixes._binary_significant_features(
            radar, mask_field, size_bins=size_bins, size_limits=size_limits,
            structure=structure, debug=debug, verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(mask_field, 1, op='and')

    return gatefilter


def echo_boundaries(
        radar, gatefilter=None, texture_window=(3, 3), texture_sample=5,
        min_texture=None, bounds_percentile=95.0, remove_small_features=False,
        size_bins=75, size_limits=(0, 300), rays_wrap_around=False,
        fill_value=None, sqi_field=None, text_field=None, bounds_field=None,
        debug=False, verbose=False):
    """
    Objectively determine the location of significant echo boundaries through
    analysis of signal quality index (SQI) texture. The a priori assumption is
    that at echo boundaries (e.g., cloud boundaries), the SQI field decreases
    substantially and therefore the SQI texture field is large near echo
    boundaries.

    Parameters
    ----------
    radar : Radar
        Radar object containing the SQI field used to derive significant echo
        boundaries.

    Returns
    -------
    gatefilter : GateFilter

    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if sqi_field is None:
        sqi_field = get_field_name('normalized_coherent_power')
    if text_field is None:
        text_field = '{}_texture'.format(sqi_field)
    if bounds_field is None:
        bounds_field = 'echo_boundaries_mask'

    if verbose:
        print 'Computing significant echo boundaries mask'

    # Compute signal quality index texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, sqi_field, ray_window=ray_window, gate_window=gate_window,
        min_sample=texture_sample, min_ncp=None, min_sweep=None,
        max_sweep=None, rays_wrap_around=rays_wrap_around,
        fill_value=fill_value, text_field=text_field, ncp_field=None)

    if min_texture is None:

        # The specified boundary percentile defines the minimum SQI texture
        # value for significant echo boundaries
        min_texture = np.percentile(
            radar.fields[text_field]['data'].compressed(), bounds_percentile,
            overwrite_input=False, interpolation='linear')

        if debug:
            max_texture = radar.fields[text_field]['data'].max()
            _range = [round(min_texture, 3), round(max_texture, 3)]
            print 'Echo boundary SQI texture range: {}'.format(_range)

        # Compute percentiles for debugging purposes
        percentiles = [5, 10, 25, 50, 75, 90, 95, 99, 100]
        textures = np.percentile(
            radar.fields[text_field]['data'].compressed(), percentiles,
            overwrite_input=False, interpolation='linear')

        if debug:
            for p, texture in zip(percentiles, textures):
                print '{}% SQI texture = {:.5f}'.format(p, texture)

    # Determine radar gates which meet minimum normalized coherent power
    # texture
    is_boundary = radar.fields[text_field]['data'] >= min_texture
    is_boundary = np.ma.filled(is_boundary, False)

    # Create significant echo boundaries field dictionary
    bounds_dict = {
        'data': is_boundary.astype(np.int8),
        'standard_name': bounds_field,
        'long_name': 'Significant echo boundaries mask',
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = not an echo boundary, 1 = echo boundary',
    }
    radar.add_field(bounds_field, bounds_dict, replace_existing=True)

    # Remove insignificant features from significant echo boundaries mask
    if remove_small_features:
        basic_fixes._binary_significant_features(
            radar, bounds_field, size_bins=size_bins, size_limits=size_limits,
            structure=structure, debug=debug, verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(bounds_field, 1, op='and')

    return gatefilter


def velocity_coherency(
        radar, gatefilter=None, text_bins=40, text_limits=(0, 20),
        nyquist=None, texture_window=(3, 3), texture_sample=5,
        max_texture=None, rays_wrap_around=False, remove_small_features=False,
        size_bins=75, size_limits=(0, 300), fill_value=None, vdop_field=None,
        text_field=None, coherent_field=None, debug=False, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if vdop_field is None:
        vdop_field = get_field_name('velocity')
    if text_field is None:
        text_field = '{}_texture'.format(vdop_field)
    if coherent_field is None:
        coherent_field = '{}_coherency_mask'.format(vdop_field)

    if verbose:
        print 'Computing Doppler velocity coherency mask'

    # Parse Nyquist velocity
    if nyquist is None:
        nyquist = radar.get_nyquist_vel(0, check_uniform=True)

    if debug:
        print 'Radar Nyquist velocity: {:.3f} m/s'.format(nyquist)

    # Compute Doppler velocity texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, vdop_field, ray_window=ray_window, gate_window=gate_window,
        min_sample=texture_sample, min_ncp=None, min_sweep=None,
        max_sweep=None, fill_value=fill_value, ncp_field=None,
        text_field=text_field)

    # Automatically determine the maximum Doppler velocity texture value which
    # brackets the coherent part of the Doppler velocity texture distribution
    if max_texture is None:

        # Bin and count Doppler velocity texture field
        counts, bin_edges = np.histogram(
            radar.fields[text_field]['data'].compressed(), bins=text_bins,
            range=text_limits, normed=False, weights=None, density=False)

        # Compute bin centers and bin width
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
        bin_width = np.diff(bin_edges).mean()

        if debug:
            print 'Bin width: {:.3f} m/s'.format(bin_width)

        # Determine the location of extrema in the Doppler velocity texture
        # distribution
        kmin = argrelextrema(counts, np.less, order=1, mode='clip')[0]
        kmax = argrelextrema(counts, np.greater, order=1, mode='clip')[0]

        if debug:
            print 'Minima located at: {} m/s'.format(bin_centers[kmin])
            print 'Maxima located at: {} m/s'.format(bin_centers[kmax])

        # Compute the theoretical noise peak location from Guassian noise
        # statistics
        noise_peak_theory = 2.0 * nyquist / np.sqrt(12.0)

        if debug:
            print 'Theoretical noise peak: {:.3f} m/s'.format(
                noise_peak_theory)

        # Find the closest Doppler velocity texture distribution peak to the
        # computed theoretical location
        # Here we assume that the Doppler velocity texture distribution
        # has at least one primary mode which correspondes to the incoherent
        # (noisy) part of the Doppler velocity texture distribution
        # Depending on the radar volume and the bin width used to define the
        # distribution, the distribution may be bimodal, with the new peak
        # corresponding to the coherent part of the Doppler velocity texture
        # distribution
        idx = np.abs(bin_centers[kmax] - noise_peak_theory).argmin()
        noise_peak = bin_centers[kmax][idx]

        if debug:
            print 'Computed noise peak: {:.3f} m/s'.format(noise_peak)

        # Determine primary and secondary peak locations for debugging
        # purposes
        if kmax.size > 1:
            counts_max = np.sort(counts[kmax], kind='mergesort')[::-1]
            prm_peak = bin_centers[np.abs(counts - counts_max[0]).argmin()]
            sec_peak = bin_centers[np.abs(counts - counts_max[1]).argmin()]

            if debug:
                print 'Primary peak: {:.3f} m/s'.format(prm_peak)
                print 'Secondary peak: {:.3f} m/s'.format(sec_peak)

        # Determine the left edge of the noise mode
        # Where this mode becomes a minimum defines the separation between
        # coherent and incoherent Doppler velocity texture values
        # TODO: if the chosen bin width is very small than multiple extrema can
        # exist such that the first minimum to the left of the noise peak is
        # not the appropriate minimum and the left edge detection breaks down
        is_left_side = bin_centers[kmin] < noise_peak
        max_texture = bin_centers[kmin][is_left_side].max() + bin_width / 2.0

        if debug:
            _range = [0.0, round(max_texture, 3)]
            print 'Doppler velocity coherent mode: {} m/s'.format(_range)

    # Create the Doppler velocity texture coherency mask
    is_coherent = np.logical_and(
            radar.fields[text_field]['data'] >= 0.0,
            radar.fields[text_field]['data'] <= max_texture)
    is_coherent = np.ma.filled(is_coherent, False)

    coherent_dict = {
        'data': is_coherent.astype(np.int8),
        'long_name': 'Doppler velocity coherency mask',
        'standard_name': coherent_field,
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = incoherent velocity, 1 = coherent velocity',
    }
    radar.add_field(coherent_field, coherent_dict, replace_existing=True)

    # Remove insignificant features from Doppler velocity coherency mask
    if remove_small_features:
        basic_fixes._binary_significant_features(
            radar, coherent_field, size_bins=size_bins,
            size_limits=size_limits, structure=structure, debug=debug,
            verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(coherent_field, 1, op='and')

    return gatefilter


def velocity_phasor_coherency(
        radar, gatefilter=None, text_bins=40, text_limits=(0, 20),
        nyquist=None, texture_window=(3, 3), texture_sample=5,
        max_texture=None, rays_wrap_around=False, remove_small_features=False,
        size_bins=75, size_limits=(0, 300), fill_value=None, vdop_field=None,
        phasor_field=None, text_field=None, coherent_field=None, debug=False,
        verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if vdop_field is None:
        vdop_field = get_field_name('velocity')
    if phasor_field is None:
        phasor_field = '{}_phasor_real'.format(vdop_field)
    if text_field is None:
        text_field = '{}_texture'.format(phasor_field)
    if coherent_field is None:
        coherent_field = '{}_coherency_mask'.format(phasor_field)

    if verbose:
        print 'Computing Doppler velocity phasor coherency mask'

    # Parse Nyquist velocity
    if nyquist is None:
        nyquist = radar.get_nyquist_vel(0, check_uniform=True)

    if debug:
        print 'Radar Nyquist velocity: {:.3f} m/s'.format(nyquist)

    # Compute the real part of Doppler velocity phasor
    # Normalize real part of phasor to the Nyquist interval
    vdop = radar.fields[vdop_field]['data']
    phasor_real = nyquist * np.cos(np.radians(360.0 * vdop / nyquist))

    # Mask invalid values
    phasor_real = np.ma.masked_invalid(phasor_real)
    phasor_real.set_fill_value(fill_value)

    # Create Doppler velocity phasor field dictionary
    phasor_dict = {
        'data': phasor_real.astype(np.float32),
        'long_name': 'Real part of Doppler velocity phasor',
        'standard_name': phasor_field,
        'valid_min': -nyquist,
        'valid_max': nyquist,
        '_FillValue': phasor_real.fill_value,
        'units': 'meters_per_second',
        'comment': ('Real part of Doppler velocity phasor normalized to the '
                    'Nyquist interval'),
    }
    radar.add_field(phasor_field, phasor_dict, replace_existing=True)

    # Compute Doppler velocity phasor texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, phasor_field, ray_window=ray_window, gate_window=gate_window,
        min_sample=texture_sample, min_ncp=None, min_sweep=None,
        max_sweep=None, fill_value=fill_value, ncp_field=None)

    # Automatically bracket coherent part of Doppler velocity phasor texture
    # distribution
    if max_texture is None:

        # Bin Doppler velocity phasor texture data and count occurrences
        # Compute bin centers and bin width
        counts, bin_edges = np.histogram(
            radar.fields[text_field]['data'].compressed(), bins=text_bins,
            range=text_limits, normed=False, weights=None, density=False)
        bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
        bin_width = np.diff(bin_edges).mean()

        if debug:
            print 'Bin width: {:.3f} m/s'.format(bin_width)

        # Determine positions of the extrema in the Doppler velocity phasor
        # texture distribution
        kmin = argrelextrema(counts, np.less, order=1, mode='clip')[0]
        kmax = argrelextrema(counts, np.greater, order=1, mode='clip')[0]

        if debug:
            print 'Minima located at: {} m/s'.format(bin_centers[kmin])
            print 'Maxima located at: {} m/s'.format(bin_centers[kmax])

        # Compute the theoretical noise peak location from Guassian noise
        # statistics
        noise_peak_theory = 2.0 * nyquist / np.sqrt(12.0)

        if debug:
            print 'Theoretical noise peak: {:.3f} m/s'.format(
                noise_peak_theory)

        # Find the closest Doppler velocity phasor texture distribution peak to
        # the computed theoretical location
        # Here we assume that the Doppler velocity phasor texture distribution
        # has at least one primary mode which correspondes to the incoherent
        # (noisy) part of the Doppler velocity phasor texture distribution
        # Depending on the radar volume and the bin width used to define the
        # distribution, the distribution may be bimodal, with the new peak
        # corresponding to the coherent part of the Doppler velocity phasor
        # texture distribution
        idx = np.abs(bin_centers[kmax] - noise_peak_theory).argmin()
        noise_peak = bin_centers[kmax][idx]

        if debug:
            print 'Computed noise peak: {:.3f} m/s'.format(noise_peak)

        # Determine primary and secondary peak locations for debugging
        # purposes
        if kmax.size > 1:
            counts_max = np.sort(counts[kmax], kind='mergesort')[::-1]
            prm_peak = bin_centers[np.abs(counts - counts_max[0]).argmin()]
            sec_peak = bin_centers[np.abs(counts - counts_max[1]).argmin()]

            if debug:
                    print 'Primary peak: {:.3f} m/s'.format(prm_peak)
                    print 'Secondary peak: {:.3f} m/s'.format(sec_peak)

        # Determine the left edge of the noise distribution
        # Where this distribution becomes a minimum will define the separation
        # between coherent and incoherent Doppler velocity phasor values
        is_left_side = bin_centers[kmin] < noise_peak
        max_texture = bin_centers[kmin][is_left_side].max() + bin_width / 2.0

        if debug:
            _range = [0.0, round(max_texture, 3)]
            print 'Doppler velocity phasor coherency mode: {} m/s'.format(
                _range)

    # Create Doppler velocity phasor coherency mask
    is_coherent = np.logical_and(
        radar.fields[text_field]['data'] >= 0.0,
        radar.fields[text_field]['data'] <= max_texture)
    is_coherent = np.ma.filled(is_coherent, False)

    # Create Doppler velocity phasor coherency mask dictionary
    coherent_dict = {
        'data': is_coherent.astype(np.int8),
        'long_name': 'Doppler velocity phasor coherency',
        'standard_name': coherent_field,
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = incoherent value, 1 = coherent value',
    }
    radar.add_field(coherent_field, coherent_dict, replace_existing=True)

    # Remove insignificant features from Doppler velocity phasor coherency mask
    if remove_small_features:
        basic_fixes._binary_significant_features(
            radar, coherent_field, size_bins=size_bins,
            size_limits=size_limits, structure=structure, debug=debug,
            verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(coherent_field, 1, op='and')

    return gatefilter


def _spectrum_width_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, rays_wrap_around=False, remove_salt=False,
        salt_window=(5, 5), salt_sample=10, fill_value=None, width_field=None,
        width_text_field=None, cohere_field=None, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if width_field is None:
        width_field = get_field_name('spectrum_width')
    if width_text_field is None:
        width_text_field = '{}_texture'.format(width_field)
    if cohere_field is None:
        cohere_field = '{}_coherency_mask'.format(width_field)

    # Compute spectrum width texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, width_field, ray_window=ray_window, gate_window=gate_window,
        min_sample=texture_sample, min_ncp=None, min_sweep=None,
        max_sweep=None, fill_value=fill_value, ncp_field=None)

    # Automatically bracket noise distribution
    if min_sigma is None and max_sigma is None:

        # Compute spectrum width texture frequency counts
        # Normalize frequency counts and compute bin centers and bin width
        width_sigma = radar.fields[width_text_field]['data']
        hist, edges = np.histogram(
            width_sigma.compressed(), bins=num_bins, range=limits,
            normed=False, weights=None, density=False)
        hist = hist.astype(np.float64) / hist.max()
        width = np.diff(edges).mean()
        half_width = width / 2.0
        bins = edges[:-1] + half_width

        if verbose:
            print 'Bin width = %.2f m/s' % width

        # Determine distribution extrema locations
        k_min = argrelextrema(
            hist, np.less, axis=0, order=1, mode='clip')[0]
        k_max = argrelextrema(
            hist, np.greater, axis=0, order=1, mode='clip')[0]

        if verbose:
            print 'Minima located at %s m/s' % bins[k_min]
            print 'Maxima located at %s m/s' % bins[k_max]

        #  Potentially a clear air volume
        if k_min.size <= 1 or k_max.size <= 1:

            # Bracket noise distribution
            # Add (left side) or subtract (right side) the half bin width to
            # account for the bin width
            max_sigma = bins.max() + half_width

            # Account for the no coherent signal case
            if k_min.size == 0:
                min_sigma = bins.min() - half_width
            else:
                min_sigma = bins[k_min][0] + half_width

            if verbose:
                print 'Computed min_sigma = %.2f m/s' % min_sigma
                print 'Computed max_sigma = %.2f m/s' % max_sigma
                print 'Radar volume is likely a clear air volume'

        # Typical volume containing sufficient scatterers (e.g., hydrometeors,
        # insects, etc.)
        else:

            # Compute primary and secondary peak locations
            hist_max = np.sort(hist[k_max], kind='mergesort')[::-1]
            prm_peak = bins[np.abs(hist - hist_max[0]).argmin()]
            sec_peak = bins[np.abs(hist - hist_max[1]).argmin()]

            if verbose:
                print 'Primary peak located at %.2f m/s' % prm_peak
                print 'Secondary peak located at %.2f m/s' % sec_peak

            # If the primary (secondary) peak velocity texture is greater than
            # the secondary (primary) peak velocity texture, than the primary
            # (secondary) peak defines the noise distribution
            noise_peak = np.max([prm_peak, sec_peak])

            if verbose:
                print 'Noise peak located at %.2f m/s' % noise_peak

            # Determine left/right sides of noise distribution
            left_side = bins[k_min] < noise_peak
            right_side = bins[k_min] > noise_peak

            # Bracket noise distribution
            # Add (left side) or subtract (right side) the half bin width to
            # account for the bin width
            min_sigma = bins[k_min][left_side].max() + half_width
            max_sigma = bins.max() + half_width

            if verbose:
                print 'Computed min_sigma = %.2f m/s' % min_sigma
                print 'Computed max_sigma = %.2f m/s' % max_sigma

    # Create the spectrum width texture coherency mask
    mask = np.logical_or(
        radar.fields[width_text_field]['data'] <= min_sigma,
        radar.fields[width_text_field]['data'] >= max_sigma)
    mask = np.ma.filled(mask, False)

    mask_dict = {
        'data': mask.astype(np.int8),
        'long_name': 'Spectrum width coherency mask',
        'standard_name': cohere_field,
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': ('0 = incoherent spectrum width, '
                    '1 = coherent spectrum width'),
    }
    radar.add_field(cohere_field, mask_dict, replace_existing=True)

    # Remove salt and pepper noise from mask
    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=[cohere_field], salt_window=salt_window,
            salt_sample=salt_sample, rays_wrap_around=rays_wrap_around,
            fill_value=0, mask_data=False, verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(cohere_field, 1, op='and')

    return gatefilter


def _significant_features(
        radar, field, gatefilter=None, min_size=None, size_bins=75,
        size_limits=(0, 300), structure=None, remove_size_field=True,
        fill_value=None, size_field=None, debug=False, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if size_field is None:
        size_field = '{}_feature_size'.format(field)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Parse binary structuring element
    if structure is None:
        structure = ndimage.generate_binary_structure(2, 1)

    # Initialize echo feature size array
    size_data = np.zeros_like(
        radar.fields[field]['data'], subok=False, dtype=np.int32)

    # Loop over all sweeps
    feature_sizes = []
    for sweep in radar.iter_slice():

        # Parse radar sweep data and define only valid gates
        is_valid_gate = ~radar.fields[field]['data'][sweep].mask

        # Label the connected features in radar sweep data and create index
        # array which defines each unique label (feature)
        labels, nlabels = ndimage.label(
            is_valid_gate, structure=structure, output=None)
        index = np.arange(1, nlabels + 1, 1)

        if debug:
            print 'Number of unique features for {}: {}'.format(sweep, nlabels)

        # Compute the size (in radar gates) of each echo feature
        # Check for case where no echo features are found, e.g., no data in
        # sweep
        if nlabels > 0:
            sweep_sizes = ndimage.labeled_comprehension(
                is_valid_gate, labels, index, np.count_nonzero, np.int32, 0)
            feature_sizes.append(sweep_sizes)

            # Set each label (feature) to its total size (in radar gates)
            for label, size in zip(index, sweep_sizes):
                size_data[sweep][labels == label] = size

    # Stack sweep echo feature sizes
    feature_sizes = np.hstack(feature_sizes)

    # Compute histogram of echo feature sizes, bin centers and bin
    # width
    counts, bin_edges = np.histogram(
        feature_sizes, bins=size_bins, range=size_limits, normed=False,
        weights=None, density=False)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
    bin_width = np.diff(bin_edges).mean()

    if debug:
        print 'Bin width: {} gate(s)'.format(bin_width)

    # Compute the peak of the echo feature size distribution
    # We expect the peak of the echo feature size distribution to be close to 1
    # radar gate
    peak_size = bin_centers[counts.argmax()] - bin_width / 2.0

    if debug:
        print 'Feature size at peak: {} gate(s)'.format(peak_size)

    # Determine the first instance when the count (sample size) for an echo
    # feature size bin reaches 0 after the distribution peak
    # This will define the minimum echo feature size
    is_zero_size = np.logical_and(
        bin_centers > peak_size, np.isclose(counts, 0, atol=1.0e-5))
    min_size = bin_centers[is_zero_size].min() - bin_width / 2.0

    if debug:
        _range = [0.0, min_size]
        print 'Insignificant feature size range: {} gates'.format(_range)

    # Mask invalid feature sizes, e.g., zero-size features
    size_data = np.ma.masked_equal(size_data, 0, copy=False)
    size_data.set_fill_value(fill_value)

    # Add echo feature size field to radar
    size_dict = {
        'data': size_data.astype(np.int32),
        'standard_name': size_field,
        'long_name': '',
        '_FillValue': size_data.fill_value,
        'units': 'unitless',
    }
    radar.add_field(size_field, size_dict, replace_existing=True)

    # Update gate filter
    gatefilter.include_above(size_field, min_size, op='and', inclusive=False)

    # Remove eacho feature size field
    if remove_size_field:
        radar.fields.pop(size_field, None)

    return gatefilter
