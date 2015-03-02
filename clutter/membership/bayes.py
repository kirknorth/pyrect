"""
clutter.membership.bayes
========================
"""

import pickle
import numpy as np

from pyart.config import get_fillvalue, get_field_name

from . import member



def _cond_prob(xi, pdf, bins, zero=1.0e-10):
    """
    Parse the conditional probability of the input data given its probability
    density function.
    """

    # Compute the probability of xi occuring given the underlying probability
    # density
    P = pdf[np.abs(bins - xi).argmin()]

    # Account for the zero frequency (probability) situation
    if P < zero:
        P = zero

    return P


def classify(radar, textures=None, moments=None, clutter_map=None, weights=1.0,
             class_prob='equal', zero=1.0e-10, fill_value=None,
             refl_field=None, precip_field=None, clutter_field=None,
             verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if refl_field is None:
        refl_field = get_field_name('reflectivity')
    if precip_field is None:
        precip_field = 'precipitation'
    if clutter_field is None:
        clutter_field = 'non-precipitation'

    # Check if at least one input is available
    if textures is None and moments is None and clutter_map is None:
        raise ValueError('No inputs specified')

    # Check if precipitation and non-precipitation fields are present for
    # textures and moments

    # TODO: parse number of inputs for each class
    n = 1.0

    # Parse class probability P(c)
    if class_prob.upper() == 'EQUAL':
        P_c = 0.5

    # Parse input probabilities P(x1, x2, ... , xn)
    if isinstance(weights, float):
        P_xi = weights

    # Initialize total probability array
    P_tot = {label: np.ones_like(radar.fields[refl_field]['data'])
             for label in textures.keys()}

    # Process radar texture fields
    for label, texture in textures.iteritems():
        for field, histogram in texture.iteritems():

            # Parse radar data and its texture distribution
            data = np.ma.filled(radar.fields[field]['data'], fill_value)
            pdf = histogram['probability density']
            bins = histogram['bin centers']

            # Compute conditional probability for each radar gate
            P_cond = member.conditional_all(
                data, pdf, bins, zero=zero, fill_value=fill_value)

            # Mask invalid probabilities
            P_cond = np.ma.masked_equal(P_cond, fill_value)

            # Bayes classifier
            P_tot[label] *= P_c * P_cond / P_xi

    # Process radar moments
    if moments is not None:
        for label, moment in moments.iteritems():
            for field, histogram in moment.iteritems():
                continue

    # Process clutter frequency map

    # Determine where the clutter (non-precipitating) class has the highest
    # probability
    is_clutter = P_tot[clutter_field] >= P_tot[precip_field]
    is_clutter = np.ma.filled(is_clutter, -1)

    mask = {
        'data': is_clutter.astype(np.int32),
        'long_name': 'Clutter classification',
        'standard_name': 'clutter_classification',
        '_FillValue': None,
        'units': None,
        'comment': ('-1 = Missing gate, 0 = Valid echo (precipitation), '
                    '1 = Clutter (non-precipitation)')
    }
    radar.add_field('clutter_classification', mask, replace_existing=True)

    return
