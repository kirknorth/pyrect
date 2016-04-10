"""
pyrect.util.texture
===================

Submodule for computing textures of radar measurements. A texture is defined
as the standard deviation of a radar measurement within a 1-D or 2-D window
centered around a radar gate.

.. autosummary::
    :toctree: generated/

    add_textures
    _compute_texture

"""

import numpy as np
from scipy import ndimage

from pyart.config import get_field_name, get_metadata, get_fillvalue


def add_textures(radar, fields=None, gatefilter=None, size=(3, 3),
                 rays_wrap_around=False, debug=False, verbose=False):
    """
    Add texture fields to radar.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar containing specified fields.
    fields : sequence of strs, optional
        Radar fields. If None, all radar fields are used.
    gatefilter : pyart.filters.GateFilter, optional
        GateFilter specifying radar gates to exclude when computing texture.
    size : int or sequence of ints, optional
        The sizes of the texture filter are given for each axis as a sequence,
        or as a single number, in which case the size is equal for all axes.
        The default filter is 3 rays and 3 gates in size.

    Optional parameters
    -------------------
    rays_wrap_around : bool, optional
        True if rays are contiguous in all sweeps.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    # parse fields
    if fields is None:
        fields = radar.fields.keys()

    for field in fields:
        if verbose:
            print 'Computing texture field: {}'.format(field)

        _compute_texture(
            radar, field, gatefilter=gatefilter, size=size,
            rays_wrap_around=rays_wrap_around, debug=debug, verbose=verbose)

    return


def _compute_texture(radar, field, gatefilter=None, size=(3, 3),
                     rays_wrap_around=False, debug=False, verbose=False):
    """
    Compute texture field of radar measurement.

    Parameters
    ----------
    radar : pyart.core.Radar
        Radar containing specified field.
    field : str
        Radar field.
    gatefilter : pyart.filters.GateFilter, optional
        GateFilter specifying radar gates to exclude when computing texture.
    size : int or sequence of ints, optional
        The sizes of the texture filter are given for each axis as a sequence,
        or as a single number, in which case the size is equal for all axes.
        The default filter is 3 rays and 3 gates in size.

    Optional parameters
    -------------------
    rays_wrap_around : bool, optional
        True if rays are contiguous in all sweeps.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Notes
    -----

    """

    # parse image data
    image = radar.fields[field]['data']

    # Parse gate filter
    if gatefilter is not None:
        image = np.ma.masked_where(gatefilter.gate_excluded, image)

    if debug:
        N = np.ma.count(image)
        print 'Sample size of data: {}'.format(N)

    # parse image border parameter
    if rays_wrap_around:
        mode = 'wrap'
    else:
        mode = 'reflect'

    x = np.empty_like(image)  # store mean of image
    y = np.empty_like(image)  # store mean of image squared
    for slc in radar.iter_slice():
        ndimage.uniform_filter(
            image[slc], size=size, output=x[slc], mode=mode, origin=0)
        ndimage.uniform_filter(
            image[slc]**2.0, size=size, output=y[slc], mode=mode, origin=0)
    std = np.ma.sqrt(y - x**2.0)

    if debug:
        N = np.ma.count(std)
        print 'Sample size of texture field: {}'.format(N)

    # parse fill value
    fill_value = radar.fields[field].get('_FillValue', get_fillvalue())
    if np.ma.is_masked(std):
        std.set_fill_value(fill_value)

    std_dict = get_metadata(field)
    std_dict['data'] = std.astype(np.float32)
    std_dict['_FillValue'] = fill_value
    radar.add_field(
        '{}_texture'.format(field), std_dict, replace_existing=True)

    return
