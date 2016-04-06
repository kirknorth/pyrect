"""
echo.correct.binary
===================

Submodule for processing binary field data.

"""

import numpy as np

from scipy import ndimage

# Necessary and/or potential future improvements to the binary submodule:
#
# * Properly handle contiguous sweep radar volumes, e.g., 360 deg PPI radar
#   volumes.


def label_size(radar, field, structure=None, rays_wrap_around=False,
               size_field=None, debug=False, verbose=False):
    """
    Label connected pixels in a binary radar field and compute the size of each
    label. Unexpected results may occur if the specified field is not binary.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing the specified field to be labeled.
    field : str
        Radar field to be labeled.
    structure : array_like, optional
        Binary structuring element used in labeling. The default structuring
        element has a squared connectivity equal to 1, i.e., only nearest
        neighbours are connected to the structure origin and
        diagonally-connected elements are not considered neighbours.
    rays_wrap_around : bool, optional
        True if all sweeps have contiguous rays, e.g., 360 deg PPI radar
        volumes.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    if verbose:
        print 'Performing binary labeling: {}'.format(field)

    if size_field is None:
        size_field = '{}_feature_size'.format(field)

    size_dict = {
        'data': np.zeros_like(radar.fields[field]['data'], dtype=np.int32),
        'units': 'unitless',
        'valid_min': 0,
        'comment': 'size in pixels of connected features',
        }
    radar.add_field(size_field, size_dict, replace_existing=True)

    for sweep, slc in enumerate(radar.iter_slice()):

        # Parse radar sweep data
        data = np.ma.getdata(radar.get_field(sweep, field))

        # Label connected pixels defined by the structuring element
        labels, nlabels = ndimage.label(data, structure=structure)
        index = np.arange(1, nlabels + 1)
        if debug:
            print 'Unique features in sweep {}: {}'.format(sweep, nlabels)

        if nlabels > 0:

            # Compute the size in number of pixels of each labeled feature
            sizes = ndimage.labeled_comprehension(
                data, labels, index, np.count_nonzero, np.int32, 0)

            # Set each labeled feature to its total size in radar gates
            for label, size in zip(index, sizes):
                radar.fields[size_field]['data'][slc][labels == label] = size

    return


def dilation(radar, field, structure=None, iterations=1,
             rays_wrap_around=False, debug=False, verbose=False):
    """
    Use binary dilation to expand around the edges of a binary radar field.
    Non-zero elements of the radar field form the subset to be dilated.
    Unexpected results may occur if the specified field is not binary.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing the specified field to be dilated.
    field : str
        Radar field to be dilated.
    structure : array_like, optional
        Binary structuring element used in dilation. The default structuring
        element has a squared connectivity equal to 1, i.e., only nearest
        neighbours are connected to the structure origin and
        diagonally-connected elements are not considered neighbours.
    iterations : int, optional
        The number of iterations to repeat dilation. If iterations is less than
        1, the dilation is repeated until the result does not change anymore.
    rays_wrap_around : bool, optional
        True if all sweeps have contiguous rays, e.g., 360 deg PPI radar
        volumes.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    if verbose:
        print 'Performing binary dilation: {}'.format(field)

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size before dilation: {}'.format(N)

    for sweep, slc in enumerate(radar.iter_slice()):

        # Parse radar sweep data
        data = np.ma.getdata(radar.get_field(sweep, field))

        # Dilate binary data
        data = ndimage.binary_dilation(
            data, structure=structure, iterations=iterations, mask=None,
            border_value=0, origin=0, brute_force=False)

        radar.fields[field]['data'][slc] = data

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size after dilation: {}'.format(N)

    return


def erosion(radar, field, structure=None, iterations=1, rays_wrap_around=False,
            debug=False, verbose=False):
    """
    Use binary erosion to shrink the edges of a binary radar field. Non-zero
    elements of the radar field form the subset to be eroded. Unexpected
    results may occur if the specified field is not binary.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing the specified field to be eroded.
    field : str
        Radar field to be eroded.
    structure : array_like, optional
        Binary structuring element used in erosion. The default structuring
        element has a squared connectivity equal to 1, i.e., only nearest
        neighbours are connected to the structure origin and
        diagonally-connected elements are not considered neighbours.
    iterations : int, optional
        The number of iterations to repeat erosion. If iterations is less than
        1, the erosion is repeated until the result does not change anymore.
    rays_wrap_around : bool, optional
        True if all sweeps have contiguous rays, e.g., 360 deg PPI radar
        volumes.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    if verbose:
        print 'Performing binary erosion: {}'.format(field)

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size before erosion: {}'.format(N)

    for sweep, slc in enumerate(radar.iter_slice()):

        # Parse radar sweep data
        data = np.ma.getdata(radar.get_field(sweep, field))

        # Erode binary data
        data = ndimage.binary_erosion(
            data, structure=structure, iterations=iterations, mask=None,
            border_value=0, origin=0, brute_force=False)

        radar.fields[field]['data'][slc] = data

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size after erosion: {}'.format(N)

    return


def fill_holes(radar, field, structure=None, rays_wrap_around=False,
               debug=False, verbose=False):
    """
    Fill holes found in a binary radar field. Non-zero elements of the radar
    field form the subset to be filled. Unexpected results may occur if the
    specified field is not binary.

    Parameters
    ----------
    radar : Radar
        Py-ART Radar containing the specified field to be eroded.
    field : str
        Radar field to be filled.
    structure : array_like, optional
        Binary structuring element used in filling. The default structuring
        element has a squared connectivity equal to 1, i.e., only nearest
        neighbours are connected to the structure origin and
        diagonally-connected elements are not considered neighbours.
    rays_wrap_around : bool, optional
        True if all sweeps have contiguous rays, e.g., 360 deg PPI radar
        volumes.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    if verbose:
        print 'Performing binary fill holes: {}'.format(field)

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size before filling holes: {}'.format(N)

    for sweep, slc in enumerate(radar.iter_slice()):

        # Parse radar sweep data
        data = np.ma.getdata(radar.get_field(sweep, field))

        # Fill binary data
        data = ndimage.binary_fill_holes(data, structure=structure, origin=0)

        radar.fields[field]['data'][slc] = data

    if debug:
        N = np.count_nonzero(radar.fields[field]['data'])
        print 'Sample size after filling holes: {}'.format(N)

    return
