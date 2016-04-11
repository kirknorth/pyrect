""" Unit tests for texture.py module. """

import numpy as np

from pyart.filters import GateFilter
from pyart.testing import sample_objects
from pyart.config import get_metadata, get_field_name

from pyrect.util import texture


def test_texture_constant_field():
    radar = _make_constant_refl_radar()
    texture._compute_texture(radar, 'reflectivity', gatefilter=None)
    refl_texture = radar.fields['reflectivity_texture']['data']

    assert np.allclose(refl_texture, 0.0, atol=1.0e-5)


def test_texture_all_excluded():
    radar = _make_constant_refl_radar()
    gatefilter = GateFilter(radar)
    gatefilter.exclude_all()
    texture._compute_texture(radar, 'reflectivity', gatefilter=gatefilter)
    refl_texture = radar.fields['reflectivity_texture']['data']

    assert np.all(refl_texture.mask)

    return


def _make_constant_refl_radar(fill=5.0):
    """ Create radar with constant reflectivity. """

    radar = sample_objects.make_empty_ppi_radar(101, 360, 1)
    refl_dict = get_metadata('reflectivity')
    refl_dict['data'] = np.full((radar.nrays, radar.ngates), fill)
    radar.add_field(get_field_name('reflectivity'), refl_dict)
    return radar
