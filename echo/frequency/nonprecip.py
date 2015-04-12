"""
clutter.frequency.nonprecip
===========================

"""

import os
import re
import json
import pickle
import numpy as np

from datetime import datetime

from pyart.io import read
from pyart.config import get_fillvalue, get_field_name

from ..correct import noise


def _pickle_map(nonprecip_map, filename, outdir=None):
    """
    Pickle a clutter frequency (probability) map.
    """

    # Parse output directory
    if outdir is not None:
        filename = os.path.join(outdir, filename)

    with open(filename, 'wb') as fid:
        pickle.dump(nonprecip_map, fid, protocol=pickle.HIGHEST_PROTOCOL)

    return


def map_date_range(start, stop, stamp, inpdir, date_str='[0-9]{12}',
                   date_fmt='%y%m%d%H%M%S', min_ncp=None, vcp_sweeps=None,
                   vcp_rays=None, exclude_fields=None, ncp_field=None,
                   debug=False, verbose=False):
    """
    Compute the clutter frequency (probability) map within the specified date
    range. The start and stop times should define a non-precipitating time
    period where (most) echoes present must be, by definition, clutter.

    Parameters
    ----------

    Optional Parameters
    -------------------

    Returns
    -------
    """

    # Parse field names
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Get all files with stamp in directory
    files = [os.path.join(inpdir, f) for f in sorted(os.listdir(inpdir))
             if stamp in f]

    if verbose:
        print 'Total number of radar files found = %i' % len(files)

    # Remove files outside date range
    time_str = [re.search(date_str, f).group() for f in files]
    times = [datetime.strptime(string, date_fmt) for string in time_str]
    files = [
        f for f, time in zip(files, times) if time >= start and time <= stop]

    if verbose:
        print 'Number of radar files after within date range = %i' % len(files)

    if vcp_sweeps is not None and vcp_rays is not None:
        nonprecip = np.zeros((vcp_sweeps, vcp_rays), dtype=np.float64)
    else:
        nonprecip = None

    # Loop over all files
    sample_size = 0
    for i, f in enumerate(files):

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Read radar data
        radar = read(f, exclude_fields=exclude_fields)

        # Check radar VCP
        if vcp_sweeps is not None and radar.nsweeps != vcp_sweeps:
            continue
        if vcp_rays is not None and radar.nrays != vcp_rays:
            continue

        # Initialize the non-precipitation map if not already done so
        if i == 0 and nonprecip is None:
            nonprecip = np.zeros((radar.nrays, radar.ngates), dtype=np.float64)

        # Increase sample size
        sample_size += 1

        # Find coherent pixels
        if min_ncp is not None:
            is_coherent = radar.fields[ncp_field]['data'] >= min_ncp
            is_coherent = np.ma.filled(is_coherent, Fales).astype(np.float64)
        else:
            is_coherent = np.zeros(nonprecip.shape, dtype=np.float64)

        # Find pixels that have a coherent signal
        nonprecip += is_coherent

    # Compute the probability a pixel (gate) has a valid echo during
    # non-precipitating events
    nonprecip_map = nonprecip / sample_size

    # Add clutter frequency map to radar object
    nonprecip = {
        'data': nonprecip_map,
        'long_name': 'Non-precipitating frequency map',
        'standard_name': 'nonprecip_map',
        'valid_min': 0.0,
        'valid_max': 1.0,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field('nonprecip_map', nonprecip, replace_existing=False)

    return {
        'non-precipitating map': nonprecip_map,
        'last radar': radar,
        'sample size': sample_size,
        'radar files': [os.path.basename(f) for f in files],
        'sweeps in VCP': vcp_sweeps,
        'rays in VCP': vcp_rays,
        'min NCP': min_ncp,
    }


def map_from_json(
        filename, inpdir=None, vcp_sweeps=None, vcp_rays=None, vcp_gates=None,
        min_ncp=None, use_filter=True, texture_window=(3, 3), texture_sample=5,
        vdop_bins=100, vdop_limits=(0, 20), sw_bins=50, sw_limits=(0, 5),
        remove_salt=True, salt_window=(5, 5), salt_sample=10,
        exclude_fields=None, ncp_field=None, debug=False, verbose=False):
    """
    Compute the non-precipitating frequency (probability) map from the files
    listed in a JSON file. The listed files should define a non-precipitating
    time period where (most) echoes present must be, by definition, not
    precipitation or cloud.

    Parameters
    ----------

    Optional Parameters
    -------------------

    Returns
    -------
    """

    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')

    # Parse files from JSON file
    with open(filename, 'r') as fid:
        files = json.load(fid)

    # Append input directory if provided
    if inpdir is not None:
        files = [os.path.join(inpdir, f) for f in files]

    if verbose:
        print 'Total number of radar files to process = %i' % len(files)

    # Parse non-precipitating frequency map
    if vcp_rays is not None and vcp_gates is not None:
        nonprecip = np.zeros((vcp_rays, vcp_gates), dtype=np.float64)
    else:
        nonprecip = None

    # Loop over all files
    sample_size = 0
    for i, f in enumerate(files):

        # Read radar data
        radar = read(f, exclude_fields=exclude_fields)

        # Check radar VCP parameters
        if vcp_sweeps is not None and radar.nsweeps != vcp_sweeps:
            continue
        if vcp_rays is not None and radar.nrays != vcp_rays:
            continue
        if vcp_gates is not None and radar.ngates != vcp_gates:
            continue

        if verbose:
            print 'Processing file %s' % os.path.basename(f)

        # Initialize the non-precipitation frequency map if it does not exist
        if i == 0 and nonprecip is None:
            vcp_sweeps = radar.nsweeps
            vcp_rays = radar.nrays
            vcp_gates = radar.ngates
            nonprecip = np.zeros((vcp_rays, vcp_gates), dtype=np.float64)

            if verbose:
                print 'VCP sweeps = {}'.format(vcp_sweeps)
                print 'VCP rays = {}'.format(vcp_rays)
                print 'VCP gates = {}'.format(vcp_gates)

        # Increase sample size
        sample_size += 1

        # Determine significant detection
        if use_filter:

            # Doppler velocity coherency
            gatefilter = noise.velocity_coherency(
                radar, gatefilter=None, num_bins=vdop_bins, limits=vdop_limits,
                texture_window=texture_window, texture_sample=texture_sample,
                min_sigma=None, max_sigma=None, nyquist=None,
                rays_wrap_around=False, remove_salt=remove_salt,
                salt_window=salt_window, salt_sample=salt_sample,
                fill_value=None, verbose=verbose)

            # Spectrum width coherency
            gatefilter = noise.spectrum_width_coherency(
                radar, gatefilter=gatefilter, num_bins=sw_bins,
                limits=sw_limits, texture_window=texture_window,
                texture_sample=texture_sample, min_sigma=None, max_sigma=None,
                rays_wrap_around=False, remove_salt=remove_salt,
                salt_window=salt_window, salt_sample=salt_sample,
                fill_value=None, verbose=verbose)

            # Significant detection
            gatefilter = noise.significant_detection(
                radar, gatefilter=gatefilter, remove_salt=remove_salt,
                salt_window=salt_window, salt_sample=salt_sample,
                min_ncp=min_ncp, detect_field=None, verbose=verbose)

            # Parse gate filter
            is_coherent = gatefilter.gate_included.astype(np.float64)

        elif min_ncp is not None:
            is_coherent = radar.fields[ncp_field]['data'] >= min_ncp
            is_coherent = np.ma.filled(is_coherent, False).astype(np.float64)

        else:
            raise ValueError('No way to determine significant detection')

        # Increase the non-precipitation map for all coherent gates (pixels)
        nonprecip += is_coherent

    # Compute the probability a gate (pixel) has a valid echo during
    # non-precipitating events
    nonprecip_map = nonprecip / sample_size

    # Add clutter frequency map to (last) radar object
    nonprecip = {
        'data': nonprecip_map,
        'long_name': 'Non-precipitating (clutter) frequency map',
        'standard_name': 'clutter_map',
        'valid_min': 0.0,
        'valid_max': 1.0,
        '_FillValue': None,
        'units': None,
    }
    radar.add_field('clutter_map', nonprecip, replace_existing=False)

    return {
        'non-precipitating map': nonprecip_map,
        'last radar': radar,
        'sample size': sample_size,
        'radar files': [os.path.basename(f) for f in files],
        'vcp_sweeps': vcp_sweeps,
        'vcp_rays': vcp_rays,
        'vcp_gates': vcp_gates,
        'min_ncp': min_ncp,
        'use_filter': use_filter,
        'texture_window': texture_window,
        'texture_sample': texture_sample,
        'remove_salt': remove_salt,
        'salt_window': salt_window,
        'salt_sample': salt_sample,
    }
