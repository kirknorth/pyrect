"""
clutter.correct.noise
=====================

"""

import numpy as np

from scipy.signal import argrelextrema

from pyart.correct import GateFilter
from pyart.config import get_fillvalue, get_field_name

from . import sweeps, basic_fixes
from ..texture import texture_fields


def significant_detection(
        radar, gatefilter=None, min_ncp=None, remove_salt=True,
        salt_window=(3, 3), salt_sample=5, rays_wrap_around=False,
        ncp_field=None, detect_field=None, verbose=False):
    """
    """

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if detect_field is None:
        detect_field = 'significant_detection'

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)

    # Coherent power criteria
    if min_ncp is not None:
        is_coherent = radar.fields[ncp_field]['data'] >= min_ncp
        gatefilter._merge(
            ~np.logical_or(gatefilter.gate_included, is_coherent), 'and', True)

    detect_dict = {
        'data': gatefilter.gate_included,
        'long_name': 'Radar significant detection',
        'standard_name': 'significant_detection',
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field(detect_field, detect_dict, replace_existing=True)

    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=[detect_field], salt_window=salt_window,
            salt_sample=salt_sample, rays_wrap_around=rays_wrap_around,
            fill_value=False, mask_data=False, verbose=verbose)

    # Update gate filter
    gatefilter._merge(~radar.fields[detect_field]['data'], 'new', True)

    return gatefilter


def velocity_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, nyquist=None, rays_wrap_around=False, remove_salt=True,
        salt_window=(3, 1), salt_sample=5, fill_value=None, vdop_field=None,
        vdop_text_field=None, verbose=False):
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

    # Find edges of noise curve
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

        # Determine distribution extrema locations
        k_min = argrelextrema(
            hist, np.less, axis=0, order=1, mode='clip')[0]
        k_max = argrelextrema(
            hist, np.greater, axis=0, order=1, mode='clip')[0]

        if verbose:
            print 'Minima located at %s m/s' % bins[k_min]
            print 'Maxima located at %s m/s' % bins[k_max]

        #  Definitive clear air volume with no velocity aliasing
        if k_min.size <= 1 or k_max.size <= 1:

            # Bracket noise distribution
            # Add (left side) or subtract (right side) the half bin width to
            # account for the bin width
            max_sigma = nyquist

            # Account for the no coherent signal case
            if k_min.size == 0:
                min_sigma = bins.min() - half_width
            else:
                min_sigma = bins[k_min][0] + half_width

            if verbose:
                print 'Computed min_sigma = %.2f m/s' % min_sigma
                print 'Computed max_sigma = %.2f m/s' % max_sigma
                print 'Radar volume is a clear air volume'

        # Typical volume containing sufficient scatterers (e.g., hydrometeors,
        # insects, etc.)
        # Velocity aliasing may or may not be present
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

            # Account for the no aliasing case
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

    mask_dict = {
        'data': mask.astype(np.bool),
        'long_name': 'Doppler velocity coherency',
        'standard_name': 'velocity_coherency',
        'comment': '0 = incoherent velocity, 1 = coherent velocity',
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field('velocity_coherency', mask_dict, replace_existing=True)

    # Remove salt and pepper noise from mask
    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=['velocity_coherency'], salt_window=salt_window,
            salt_sample=salt_sample, rays_wrap_around=rays_wrap_around,
            fill_value=0.0, mask_data=False, verbose=verbose)

        # Parse mask
        mask = radar.fields['velocity_coherency']['data']

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)
        gatefilter._merge(~mask, 'new', True)
    else:
        gatefilter._merge(
            ~np.logical_or(gatefilter.gate_included, mask), 'and', True)

    return gatefilter


def spectrum_width_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, rays_wrap_around=False, remove_salt=True,
        salt_window=(3, 1), salt_sample=5, fill_value=None, width_field=None,
        width_text_field=None, verbose=False):
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

    mask_dict = {
        'data': mask.astype(np.bool),
        'long_name': 'Spectrum width coherency',
        'standard_name': 'spectrum_width_coherency',
        'comment': ('0 = incoherent spectrum width, '
                    '1 = coherent spectrum width'),
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field(
        'spectrum_width_coherency', mask_dict, replace_existing=True)

    # Remove salt and pepper noise from mask
    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=['spectrum_width_coherency'],
            salt_window=salt_window, salt_sample=salt_sample,
            rays_wrap_around=rays_wrap_around, fill_value=0.0,
            mask_data=False, verbose=verbose)

        # Parse mask
        mask = radar.fields['spectrum_width_coherency']['data']

    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)
        gatefilter._merge(~mask, 'new', True)
    else:
        gatefilter._merge(
            ~np.logical_or(gatefilter.gate_included, mask), 'and', True)

    return gatefilter


def velocity_phasor_coherency(
        radar, gatefilter=None, num_bins=10, limits=None,
        texture_window=(3, 3), texture_sample=5, min_sigma=None,
        max_sigma=None, rays_wrap_around=False, remove_salt=True,
        salt_window=(3, 1), salt_sample=5, fill_value=None, vdop_field=None,
        vdop_phase_field=None, phase_text_field=None, verbose=False):
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

    # Compute the real part of Doppler velocity phasor
    vdop = radar.fields[vdop_field]['data']
    vdop_phase = np.real(np.exp(1j * np.radians(vdop)))

    # Add Doppler velocity phasor field to radar object
    phasor = {
        'data': vdop_phase,
        'long_name': '{} phasor'.format(
            radar.fields[vdop_field]['long_name']),
        'standard_name': '{}_phasor'.format(
            radar.fields[vdop_field]['standard_name']),
        '_FillValue': vdop_phase.fill_value,
        'units': None,
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

    mask_dict = {
        'data': mask.astype(np.bool),
        'long_name': 'Doppler velocity phasor coherency',
        'standard_name': 'velocity_phasor_coherency',
        'comment': ('0 = incoherent velocity phasor, '
                    '1 = coherent velocity phasor'),
        'valid_min': 0,
        'valid_max': 1,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field(
        'velocity_phasor_coherency', mask_dict, replace_existing=True)

    # Remove salt and pepper noise from mask
    if remove_salt:
        basic_fixes.remove_salt(
            radar, fields=['velocity_phasor_coherency'],
            salt_window=salt_window, salt_sample=salt_sample,
            rays_wrap_around=rays_wrap_around, fill_value=0.0,
            mask_data=False, verbose=verbose)

        # Parse mask
        mask = radar.fields['velocity_phasor_coherency']['data']


    # Parse gate filter
    if gatefilter is None:
        gatefilter = GateFilter(radar, exclude_based=False)
        gatefilter._merge(~mask, 'new', True)
    else:
        gatefilter._merge(
            ~np.logical_or(gatefilter.gate_included, mask), 'and', True)

    return gatefilter
