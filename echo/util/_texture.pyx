"""
echo.util._texture
==================

Cython routines for computing texture fields. A texture field is defined as the
standard deviation of a radar measurement within a 1-D or 2-D window centered
around a radar gate.

"""

import time
import numpy as np

cimport numpy as np
cimport cython

from cpython cimport bool

# Necessary and/or potential future improvements to _texture submodule:
#
# * Properly handle contiguous sweep volumes, e.g., 360 deg PPI volumes.

@cython.boundscheck(False)

def compute_texture(np.ndarray[np.float32_t, ndim=2] field,
                    np.ndarray[np.int32_t, ndim=1] sweep_start,
                    np.ndarray[np.int32_t, ndim=1] sweep_end,
                    np.int32_t ray_window, np.int32_t gate_window,
                    bool rays_wrap_around, np.float32_t fill_value,
                    bool debug, bool verbose):
    """

    Parameters
    ----------
    field : ndarray

    sweep_start : ndarray

    sweep_end : ndarray

    fill_value : float
        Value indicating missing or bad data in field.
    debug : bool
        True to print debugging information, False to suppress.
    verbose : bool
        True to print relevant information, False to suppress.

    Returns
    -------
    sigma : ndarray
        Texture field computed from input radar field.
    sample_size : ndarray
        Number of available radar measurements within texture window to compute
        texture field.

    """

    cdef np.int32_t nrays = field.shape[0]
    cdef np.int32_t ngates = field.shape[1]
    cdef np.int32_t nsweeps = sweep_start.shape[0]

    cdef np.ndarray[np.float32_t, ndim=2] sigma
    cdef np.ndarray[np.int32_t, ndim=2] sample_size
    cdef np.ndarray[np.int32_t, ndim=2] valid_gate

    cdef np.int32_t sweep, ray, gate, ray_start, ray_end, sample
    cdef np.int32_t r0, rf, g0, gf

    # Initialize arrays
    sigma = np.full_like(field, np.nan, dtype=np.float32)
    sample_size = np.zeros_like(field, dtype=np.int32)
    valid_gate = np.zeros_like(field, dtype=np.int32)

    #
    field[np.isclose(field, fill_value, atol=1.0e-5)] = np.nan
    valid_gate[np.isfinite(field)] = 1

    if debug:
        print('Ray window half length: {}'.format(ray_window // 2))
        print('Gate window half length: {}'.format(gate_window // 2))

    for sweep in range(nsweeps):

        # Parse the sweep start and end indices
        ray_start = sweep_start[sweep]
        ray_end = sweep_end[sweep]

        for ray in range(ray_start, ray_end + 1):

            r0 = ray - ray_window // 2
            rf = ray + ray_window // 2

            # Check for condition where current ray is close to a sweep
            # boundary
            if r0 < ray_start or rf > ray_end:

                if rays_wrap_around:
                    raise ValueError
                else:
                    if r0 < ray_start:
                        r0 = ray_start
                    if rf > ray_end:
                        rf = ray_end

            for gate in range(ngates):

                # Gates along a ray are by definition not contiguous so no
                # wrapping conditions have to be checked for
                g0 = max(0, gate - gate_window // 2)
                gf = min(ngates - 1, gate + gate_window // 2)

                # Compute sample size within texture window
                sample = np.sum(valid_gate[r0:rf+1,g0:gf+1])
                sample_size[ray,gate] = sample

                if sample > 1 and valid_gate[ray,gate]:
                    sigma[ray,gate] = np.nanstd(field[r0:rf+1,g0:gf+1], ddof=1)

    return sigma, sample_size
