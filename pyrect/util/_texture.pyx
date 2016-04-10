"""
pyrect.util._texture
====================

Cython routines for computing texture fields. A texture field is defined as the
standard deviation of a radar measurement within a 1-D or 2-D window centered
around a radar gate.

.. autosummary::
    :toctree: generated/

    add_texture

"""

import numpy as np

cimport numpy as np
cimport cython

from pyart.config import get_metadata, get_fillvalue

# Necessary and/or potential future improvements to _texture submodule:
#
# * Properly handle contiguous sweep volumes, e.g., 360 deg PPI volumes.


@cython.boundscheck(False)
@cython.wraparound(False)
def add_texture(radar, field, size=(3, 3), rays_wrap_around=False, debug=False,
                verbose=False):
    """

    Parameters
    ----------
    radar : Radar
        Radar containing specified field.
    field : str
        Radar field.
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

    cdef int sweep, ray, gate

    cdef float [:, ::1] sigma = np.empty_like(radar.fields[field]['data'])
    cdef int [:, ::1] sample
    cdef int [:, ::1] mask

    collect = _IndexCollector(radar, rays_wrap_around=rays_wrap_around)

    for ray in range(collect.nrays):
        for gate in range(collect.ngates):
            sigma[ray, gate] = 0.0

    return



cdef class _IndexCollector:
    """
    """

    cdef int nrays, ngates, nsweeps
    cdef int [:] sweep_start, sweep_end

    def __init__(self, radar, size=(3, 3), rays_wrap_around=False):
        """ Initialize. """

        # Default parameters
        self.size = size
        self.rays_wrap_around = rays_wrap_around

        self.nrays = radar.nrays
        self.ngates = radar.ngates
        self.nsweeps = radar.nsweeps
        self.sweep_start = radar.sweep_start_ray_index['data']
        self.sweep_end = radar.sweep_end_ray_index['data']

    def add_neighbors(self):
        """
        """
        return



