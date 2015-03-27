"""
echo.location.common
====================

"""


import numpy as np


def standard_refraction(radar, use_km=False):
    """
    """

    # Effective radius of Earth in meters
    Re = 6371.0 * 4.0 / 3.0 * 1000.0

    # Parse radar gate coordinates
    rng = radar.range['data']
    azi = np.radians(radar.azimuth['data'])
    ele = np.radians(radar.elevation['data'])

    # Mesh radar gate coordinates
    AZI, RNG = np.meshgrid(azi, rng, indexing='ij')
    ELE, RNG = np.meshgrid(ele, rng, indexing='ij')

    # Compute vertical height (z), arc length (s), eastward distance (x), and
    # northward distance (y) relative to radar
    z = np.sqrt(RNG**2 + 2.0 * RNG * Re * np.sin(ELE) + Re**2) - Re
    s = Re * np.arcsin(RNG * np.cos(ELE) / (z + Re))
    x = s * np.sin(AZI)
    y = s * np.cos(AZI)

    if use_km:
        x, y, z = x / 1000.0, y / 1000.0, z / 1000.0

    return x, y, z
