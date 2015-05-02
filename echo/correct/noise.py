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
        radar, gatefilter=None, remove_salt=True, salt_window=(5, 5),
        salt_sample=10, fill_holes=False, dilate=True, structure=None,
        dilate_iter=1, rays_wrap_around=False, min_ncp=None,
        ncp_field=None, detect_field=None, verbose=False):
    """
    """

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if detect_field is None:
        detect_field = 'significant_detection_mask'

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Coherent power criteria
    if min_ncp is not None and ncp_field in radar.fields:
        gatefilter.include_above(ncp_field, min_ncp, op='and')

    detect_dict = {
        'data': gatefilter.gate_included.astype(np.int8),
        'long_name': 'Radar significant detection mask',
        'standard_name': 'significant_detection_mask',
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = not significant, 1 = significant',
    }
    radar.add_field(detect_field, detect_dict, replace_existing=True)

    # Remove salt and pepper noise from significant detection mask
    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=[detect_field], salt_window=salt_window,
            salt_sample=salt_sample, rays_wrap_around=rays_wrap_around,
            fill_value=0, mask_data=False, verbose=verbose)

    # Fill holes in significant detection mask
    if fill_holes:
        basic_fixes._binary_fill(radar, detect_field, structure=structure)

    # Dilate significant detection mask
    if dilate:
        basic_fixes._binary_dilation(
            radar, detect_field, structure=structure, iterations=dilate_iter)

    # Update gate filter
    gatefilter.include_equal(detect_field, 1, op='new')

    return gatefilter


def hildebrand_noise(
        radar, gatefilter=None, scale=1.0, remove_salt=True,
        salt_window=(5, 5), salt_sample=10, rays_wrap_around=False,
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

    # Compute noise floor mask and add field to radar
    power = radar.fields[power_field]['data']
    noise = radar.fields[noise_field]['data']
    is_noise = np.ma.filled(power >= noise, False)

    mask_dict = {
        'data': is_noise.astype(np.int8),
        'long_name': 'Noise floor mask',
        'standard_name': mask_field,
        'valid_min': 0,
        'valid_max': 1,
        'units': 'unitless',
        '_FillValue': None,
        'comment': '0 = below noise floor, 1 = above noise floor',
    }
    radar.add_field(mask_field, mask_dict, replace_existing=True)

    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=[mask_field], salt_window=salt_window,
            salt_sample=salt_sample, rays_wrap_around=rays_wrap_around,
            fill_value=0, mask_data=False, verbose=verbose)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Update gate filter
    gatefilter.include_equal(mask_field, 1, op='and')

    return gatefilter


def velocity_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, nyquist=None, rays_wrap_around=False,
        remove_salt=False, salt_window=(5, 5), salt_sample=10, fill_value=None,
        vdop_field=None, vdop_text_field=None, cohere_field=None,
        verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if vdop_field is None:
        vdop_field = get_field_name('velocity')
    if vdop_text_field is None:
        vdop_text_field = '{}_texture'.format(vdop_field)
    if cohere_field is None:
        cohere_field = '{}_coherency_mask'.format(vdop_field)

    # Parse Nyquist velocity
    if nyquist is None:
        nyquist = radar.instrument_parameters['nyquist_velocity']['data'][0]

    if verbose:
        print 'Nyquist velocity = %.3f m/s' % nyquist

    # Compute Doppler velocity texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, vdop_field, ray_window=ray_window, gate_window=gate_window,
        min_sample=texture_sample, min_ncp=None, min_sweep=None,
        max_sweep=None, fill_value=fill_value, ncp_field=None)

    # Find edges of noise distribution
    if min_sigma is None and max_sigma is None:

        # Compute Doppler velocity texture frequency counts
        # Normalize frequency counts and compute bin centers and bin width
        vdop_sigma = radar.fields[vdop_text_field]['data']
        hist, edges = np.histogram(
            vdop_sigma.compressed(), bins=num_bins, range=limits, normed=False,
            weights=None, density=False)
        hist = hist.astype(np.float64) / hist.max()
        width = np.diff(edges).mean()
        half_width = width / 2.0
        bins = edges[:-1] + half_width

        if verbose:
            print 'Bin width = %.2f m/s' % width

        # Determine Doppler velocity texture distribution extrema locations
        k_min = argrelextrema(
            hist, np.less, axis=0, order=1, mode='clip')[0]
        k_max = argrelextrema(
            hist, np.greater, axis=0, order=1, mode='clip')[0]

        if verbose:
            print 'Minima located at %s m/s' % bins[k_min]
            print 'Maxima located at %s m/s' % bins[k_max]

        #  Case 1: Definitive clear air volume with no velocity aliasing
        if k_min.size <= 1 or k_max.size <= 1:

            # Bracket noise distribution
            # Add (to the left side) or subtract (to the right side) the half
            # bin width to account for the bin width
            max_sigma = nyquist

            # Case 1.1: No coherent signal in the Doppler velocity texture
            # distribution
            if k_min.size == 0:
                min_sigma = bins.min() - half_width
            else:
                min_sigma = bins[k_min][0] + half_width

            if verbose:
                print 'Computed min_sigma = %.2f m/s' % min_sigma
                print 'Computed max_sigma = %.2f m/s' % max_sigma
                print 'Radar volume is a clear air volume'

        # Case 2: Typical radar volume containing sufficient scatterers (e.g.,
        # hydrometeors, insects, etc.)
        # Velocity aliasing may or may not be present
        else:

            # Compute primary and secondary peak locations
            # One of these peaks defines the noise distribution
            hist_max = np.sort(hist[k_max], kind='mergesort')[::-1]
            prm_peak = bins[np.abs(hist - hist_max[0]).argmin()]
            sec_peak = bins[np.abs(hist - hist_max[1]).argmin()]
            peaks = np.array([prm_peak, sec_peak], dtype=np.float64)

            if verbose:
                print 'Primary peak located at %.2f m/s' % prm_peak
                print 'Secondary peak located at %.2f m/s' % sec_peak

            # The noise distribution will typically be defined as the largest
            # (in terms of Doppler velocity) of the primary and secondary
            # peaks, e.g., if the primary (secondary) peak velocity texture is
            # greater than the secondary (primary) peak velocity texture, than
            # the primary (secondary) peak defines the noise distribution
            # However in rare circumstances where the velocity aliasing signal
            # is very strong, this will not be the case
            # The most general way to determe which peak is the noise peak is
            # to use basic Guassian noise theory, which predicts that the noise
            # peak (mean) is a function of Nyquist velocity
            idx = np.abs(peaks - 2.0 * nyquist / np.sqrt(12.0)).argmin()
            noise_peak_theory = peaks[idx]
            noise_peak = np.max(peaks)

            # Case 2.1: Strong velocity aliasing signal
            if noise_peak != noise_peak_theory:
                noise_peak = noise_peak_theory
                if verbose:
                    print 'Strong velocity aliasing signal'

            if verbose:
                print 'Noise peak located at %.2f m/s' % noise_peak

            # Determine left and right sides of noise distribution
            left_side = bins[k_min] < noise_peak
            right_side = bins[k_min] > noise_peak

            # Bracket noise distribution
            # Add (to the left side) or subtract (to the right side) the half
            # bin width to account for the bin width
            min_sigma = bins[k_min][left_side].max() + half_width

            # Case 2.2: No velocity aliasing signal
            if any(right_side):
                max_sigma = bins[k_min][right_side].min() - half_width
            else:
                max_sigma = nyquist

            if verbose:
                print 'Computed min_sigma = %.2f m/s' % min_sigma
                print 'Computed max_sigma = %.2f m/s' % max_sigma

    # Create the Doppler velocity texture coherency mask
    mask = np.logical_or(
        radar.fields[vdop_text_field]['data'] <= min_sigma,
        radar.fields[vdop_text_field]['data'] >= max_sigma)
    mask = np.ma.filled(mask, False)

    mask_dict = {
        'data': mask.astype(np.int8),
        'long_name': 'Doppler velocity coherency mask',
        'standard_name': cohere_field,
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': '0 = incoherent velocity, 1 = coherent velocity',
    }
    radar.add_field(cohere_field, mask_dict, replace_existing=True)

    # Remove salt and pepper noise from Doppler velocity coherency mask
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


def spectrum_width_coherency(
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


def velocity_phasor_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, rays_wrap_around=False, remove_salt=False,
        salt_window=(5, 5), salt_sample=10, fill_value=None, vdop_field=None,
        vdop_phase_field=None, phase_text_field=None, cohere_field=None,
        verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if vdop_field is None:
        vdop_field = get_field_name('velocity')
    if vdop_phase_field is None:
        vdop_phase_field = '{}_phasor'.format(vdop_field)
    if phase_text_field is None:
        phase_text_field = '{}_texture'.format(vdop_phase_field)
    if cohere_field is None:
        cohere_field = '{}_coherency_mask'.format(vdop_phase_field)

    # Compute the real part of Doppler velocity phasor
    vdop = radar.fields[vdop_field]['data']
    vdop_phase = np.real(np.exp(1j * np.radians(vdop)))
    vdop_phase.set_fill_value(fill_value)

    # Add Doppler velocity phasor field to radar object
    phasor = {
        'data': vdop_phase.astype(np.float32),
        'long_name': 'Doppler velocity phasor',
        'standard_name': vdop_phase_field,
        'valid_min': -1.0,
        'valid_max': 1.0,
        '_FillValue': vdop_phase.fill_value,
        'units': 'unitless',
        'comment': 'Real part of phasor',
    }
    radar.add_field(vdop_phase_field, phasor, replace_existing=True)

    # Compute Doppler velocity phasor texture field
    ray_window, gate_window = texture_window
    texture_fields._compute_field(
        radar, vdop_phase_field, ray_window=ray_window,
        gate_window=gate_window, min_sample=texture_sample, min_ncp=None,
        min_sweep=None, max_sweep=None, fill_value=fill_value, ncp_field=None)

    # Automatically bracket noise distribution
    if min_sigma is None and max_sigma is None:

        # Compute Doppler velocity phasor texture frequency counts
        # Normalize histogram counts and compute bin centers and bin width
        phase_sigma = radar.fields[phase_text_field]['data']
        hist, edges = np.histogram(
            phase_sigma.compressed(), bins=num_bins, range=limits,
            normed=False, weights=None, density=False)
        hist = hist.astype(np.float64) / hist.max()
        width = np.diff(edges).mean()
        half_width = width / 2.0
        bins = edges[:-1] + half_width

        if verbose:
            print 'Bin width = %.2f' % width

        # Determine distribution extrema locations
        k_min = argrelextrema(
            hist, np.less, axis=0, order=1, mode='clip')[0]
        k_max = argrelextrema(
            hist, np.greater, axis=0, order=1, mode='clip')[0]

        if verbose:
            print 'Minima located at %s' % bins[k_min]
            print 'Maxima located at %s' % bins[k_max]

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
                print 'Computed min_sigma = %.2f' % min_sigma
                print 'Computed max_sigma = %.2f' % max_sigma
                print 'Radar volume is likely a clear air volume'

        # Typical volume containing sufficient scatterers (e.g., hydrometeors,
        # insects, etc.)
        else:

            # Compute primary and secondary peak locations
            hist_max = np.sort(hist[k_max], kind='mergesort')[::-1]
            prm_peak = bins[np.abs(hist - hist_max[0]).argmin()]
            sec_peak = bins[np.abs(hist - hist_max[1]).argmin()]

            if verbose:
                print 'Primary peak located at %.2f' % prm_peak
                print 'Secondary peak located at %.2f' % sec_peak

            # If the primary (secondary) peak velocity texture is greater than
            # the secondary (primary) peak velocity texture, than the primary
            # (secondary) peak defines the noise distribution
            noise_peak = np.max([prm_peak, sec_peak])

            if verbose:
                print 'Noise peak located at %.2f' % noise_peak

            # Determine left/right sides of noise distribution
            left_side = bins[k_min] < noise_peak
            right_side = bins[k_min] > noise_peak

            # Bracket noise distribution
            # Add (left side) or subtract (right side) the half bin width to
            # account for the bin width
            min_sigma = bins[k_min][left_side].max() + half_width
            max_sigma = bins.max() + half_width

            if verbose:
                print 'Computed min_sigma = %.2f' % min_sigma
                print 'Computed max_sigma = %.2f' % max_sigma

    # Create Doppler velocity phasor texture coherency mask
    mask = np.logical_or(
        radar.fields[phase_text_field]['data'] <= min_sigma,
        radar.fields[phase_text_field]['data'] >= max_sigma)
    mask = np.ma.filled(mask, False)

    mask_dict = {
        'data': mask.astype(np.int8),
        'long_name': 'Doppler velocity phasor coherency',
        'standard_name': cohere_field,
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': 'unitless',
        'comment': ('0 = incoherent velocity phasor, '
                    '1 = coherent velocity phasor'),
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
        radar, field, gatefilter=None, num_bins=100, limits=(0, 100),
        structure=None, fill_value=None, size_field=None, debug=False):
    """
    """

    # Parse field names
    if size_field is None:
        size_field = '{}_feature_size'.format(field)

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Parse binary structuring element
    if structure is None:
        structure = ndimage.generate_binary_structure(2, 1)

    # Initialize echo feature size field
    radar.fields[size_field] = get_metadata(size_field)
    radar.fields[size_field]['data'] = np.zeros_like(
        radar.fields[field]['data'])

    # Loop over all sweeps
    feature_sizes = []
    for sweep in radar.iter_slice():

        # Parse radar sweep data and define only valid gates
        is_valid_gate = ~radar.fields[field]['data'][sweep].mask

        # Label radar sweep data and create index array which defines each
        # unique label
        labels, nlabels = ndimage.label(
            is_valid_gate, structure=structure, output=None)
        index = np.arange(1, nlabels + 1, 1)

        if debug:
            print 'Number of unique features for {}: {}'.format(sweep, nlabels)

        # Compute the size (in radar gates) of each echo feature
        sweep_sizes = ndimage.labeled_comprehension(
            is_valid_gate, labels, index, np.sum, np.int32, 0)

        # Append sweep echo feature sizes
        feature_sizes.append(sweep_sizes)

        # Set each label to its total size (in radar gates)
        for label, size in zip(index, sweep_sizes):
            is_label = labels == label
            radar.fields[size_field]['data'][sweep][is_label] = size

    # Stack sweep echo feature sizes
    feature_sizes = np.hstack(feature_sizes)

    # Compute histogram of echo feature sizes, bin centers and bin
    # width
    counts, bin_edges = np.histogram(
        feature_sizes, bins=num_bins, range=limits, normed=False,
        weights=None, density=False)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2.0
    bin_width = np.diff(bin_edges).mean()

    if debug:
        print 'Bin width: {} gate(s)'.format(bin_width)

    # Compute the peak of the echo feature size distribution
    # We expect the peak of the echo feature size distribution to be close to 1
    # radar gate
    idx = np.abs(counts - counts.max()).argmin()
    peak_size = bin_centers[idx] - bin_width / 2.0

    if debug:
        print 'Echo feature size at peak: {} gate(s)'.format(peak_size)

    # Determine the first instance when the echo feature size reaches 0 radar
    # gates after the distribution peak
    # This will define the minimum echo feature size
    is_zero_size = np.logical_and(bin_centers > peak_size, counts == 0)
    min_size = bin_centers[is_zero_size].min() - bin_width / 2.0

    if debug:
        print 'Minimum echo feature size: {} gates'.format(min_size)

    # Update gate filter
    gatefilter.include_above(size_field, min_size, op='and', inclusive=False)

    return gatefilter
