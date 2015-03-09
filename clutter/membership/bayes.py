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
             class_prob='equal', min_inputs=1, min_ncp=0.3, zero=1.0e-10,
             ignore_inputs=None, fill_value=None, ncp_field=None,
             precip_field=None, clutter_field=None, verbose=False):
    """
    """

    # Parse fill value
    if fill_value is None:
        fill_value = get_fillvalue()

    # Parse field names
    if ncp_field is None:
        ncp_field = get_field_name('normalized_coherent_power')
    if precip_field is None:
        precip_field = 'precipitation'
    if clutter_field is None:
        clutter_field = 'non-precipitation'

    # Parse ignore fields
    if ignore_inputs is None:
        ignore_inputs = []

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

    # Initialize number of inputs array
    num_inputs = np.zeros_like(
        radar.fields[ncp_field]['data'].data, dtype=np.int32)

    # Initialize total probability and number of inputs arrays
    P_tot = {
        label: np.ones(radar.fields[ncp_field]['data'].shape, dtype=np.float64)
        for label in [precip_field, clutter_field]
        }
    num_inputs = {
        label: np.zeros(radar.fields[ncp_field]['data'].shape, dtype=np.int32)
        for label in [precip_field, clutter_field]
        }

    # Parse normalized coherent power
    ncp = radar.fields[ncp_field]['data']

    # Process radar texture fields
    if textures is not None:
        for label, texture in textures.iteritems():
            for field, histogram in texture.iteritems():

                if field in ignore_inputs:
                    continue

                # Parse radar texture data and its distribution
                data = radar.fields[field]['data']
                pdf = histogram['probability density']
                bins = histogram['bin centers']

                # Mask incoherent gates
                data = np.ma.masked_where(ncp < min_ncp, data)
                data = np.ma.filled(data, fill_value)

                # Compute conditional probability for each radar gate
                P_cond = member.conditional_all(
                data, pdf, bins, zero=zero, fill_value=fill_value)

                # Determine where conditional probability is valid
                valid_prob = P_cond != fill_value
                num_inputs[label] += valid_prob

                # Bayes classifier
                P_tot[label][valid_prob] *= P_c * P_cond[valid_prob] / P_xi

    # Process radar moments
    if moments is not None:
        for label, moment in moments.iteritems():
            for field, histogram in moment.iteritems():

                if field in ignore_inputs:
                    continue

                # Parse radar moment data and its distribution
                data = radar.fields[field]['data']
                pdf = histogram['probability density']
                bins = histogram['bin centers']

                # Mask incoherent gates
                data = np.ma.masked_where(ncp < min_ncp, data)
                data = np.ma.filled(data, fill_value)

                # Compute conditional probability for each radar gate
                P_cond = member.conditional_all(
                    data, pdf, bins, zero=zero, fill_value=fill_value)

                # Determine where conditional probability is valid
                valid_prob = P_cond != fill_value
                num_inputs[label] += valid_prob

                # Bayes classifier
                P_tot[label][valid_prob] *= P_c * P_cond[valid_prob] / P_xi

    # Process clutter frequency map
    if clutter_map is not None:
        pdf = clutter_map['clutter frequency map']

    # Mask gates where not enough inputs were available to properly classify
    for label, sample_size in num_inputs.iteritems():
        P_tot[label] = np.ma.masked_where(
            sample_size < min_inputs, P_tot[label])

    # Determine where the clutter (non-precipitating) class has the highest
    # probability
    is_clutter = P_tot[clutter_field] >= P_tot[precip_field]
    mask = np.ma.filled(is_clutter.astype(np.int32), -1)

    mask = {
        'data': mask.astype(np.int32),
        'long_name': 'Clutter classification',
        'standard_name': 'clutter_classification',
        '_FillValue': None,
        'units': None,
        'comment': ('-1 = Missing gate, 0 = Valid echo (precipitation), '
                    '1 = Clutter (non-precipitation)')
    }
    radar.add_field('clutter_classification', mask, replace_existing=True)

    return
