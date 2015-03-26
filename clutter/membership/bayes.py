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


def classify(radar, textures=None, moments=None, clutter_map=None,
             gatefilter=None, weights=1.0, class_prob='equal', min_inputs=1,
             min_ncp=None, zero=1.0e-10, ignore_inputs=None, use_insects=True,
             fill_value=None, ncp_field=None, precip_field=None,
             ground_field=None, insect_field=None, verbose=False):
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
    if ground_field is None:
        ground_field = 'ground_clutter'
    if insect_field is None:
        insect_field = 'insects'

    # Parse ignore fields
    if ignore_inputs is None:
        ignore_inputs = []

    # Check if at least one input is available
    if textures is None and moments is None and clutter_map is None:
        raise ValueError('No inputs specified')

    # Parse classification labels
    if use_insects:
        labels = [precip_field, ground_field, insect_field]
    else:
        labels = [precip_field, ground_field]
        if textures is not None:
            textures.pop(insect_field, None)
        if moments is not None:
            moments.pop(insect_field, None)

    # Determine total number of inputs available for each class
    inputs = {label: 0 for label in labels}
    if textures is not None:
        for label, texture in textures.iteritems():
            fields = [field for field in texture if field not in ignore_inputs]
            inputs[label] += len(fields)
    if moments is not None:
        for label, moment in moments.iteritems():
            fields = [field for field in moment if field not in ignore_inputs]
            inputs[label] += len(fields)
    if clutter_map is not None:
        inputs[ground_field] += 1

    if verbose:
        for label in labels:
            print 'Total number of inputs for {} = {}'.format(
                label, inputs[label])

    # Parse class probability P(c)
    if class_prob.upper() == 'EQUAL':
        P_c = 0.5

    # Parse input probabilities P(x1, x2, ... , xn)
    if isinstance(weights, float):
        P_xi = weights

    # Initialize total probability and number of inputs arrays
    P_tot = {
        label: np.ones(radar.fields[ncp_field]['data'].shape, dtype=np.float64)
        for label in labels
        }
    num_inputs = {
        label: np.zeros(radar.fields[ncp_field]['data'].shape, dtype=np.int32)
        for label in labels
        }

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
                if min_ncp is not None:
                    ncp = radar.fields[ncp_field]['data']
                    data = np.ma.masked_where(ncp < min_ncp, data)

                # Prepare data for ingest into Fortran wrapper
                data = np.ma.filled(data, fill_value)
                data = np.asfortranarray(data, dtype=np.float64)

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
                if min_ncp is not None:
                    ncp = radar.fields[ncp_field]['data']
                    data = np.ma.masked_where(ncp < min_ncp, data)

                # Prepare data for ingest into Fortran wrapper
                data = np.ma.filled(data, fill_value)
                data = np.asfortranarray(data, dtype=np.float64)

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

    # Mask excluded gates from gate filter
    if gatefilter is not None:
        for label in P_tot:
            P_tot[label] = np.ma.masked_where(
                gatefilter.gate_excluded, P_tot[label])

    # Determine where each class is most probable
    echo = np.zeros(P_tot[precip_field].shape, dtype=np.int32)
    if use_insects:
        is_ground = np.logical_and(
            P_tot[ground_field] > P_tot[precip_field],
            P_tot[ground_field] > P_tot[insect_field])
        is_insect = np.logical_and(
            P_tot[insect_field] > P_tot[ground_field],
            P_tot[insect_field] > P_tot[precip_field])
        is_missing = P_tot[ground_field].mask
        echo[is_ground] = 1
        echo[is_insect] = 2
        echo[is_missing] = -1

    else:
        is_ground = P_tot[ground_field] > P_tot[precip_field]
        is_missing = P_tot[ground_field].mask
        echo[is_ground] = 1
        echo[is_missing] = -1

    # Create echo classification dictionary and add it to the radar object
    echo = {
        'data': echo,
        'long_name': 'Echo classification',
        'standard_name': 'echo_classification',
        '_FillValue': None,
        'units': None,
        'comment': ('-1 = Missing gate, 0 = cloud or precipitation, '
                    '1 = Ground clutter, 2 = Insects')
    }
    radar.add_field('echo_classification', echo, replace_existing=True)

    return
