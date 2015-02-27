"""
clutter.membership.bayes
========================
"""

import pickle
import numpy as np

from pyart.config import get_field_name



def _conditional_prob(xi, dist, bins, zero=1.0e-10):
    """
    """

    # Compute the probability of xi from the distribution
    prob = dist(np.abs(bins - xi).argmin())

    # Account for the zero frequency (probability) situation
    if prob < zero:
        prob = zero

    return prob


def classify(radar, textures, moments=None, type='sum', weights='constant',
             zero=1.0e-10):
    """
    """

    echo = {}

    # Process radar texture fields
    for label, texture in textures.iteritems():
        for field, histogram in texture.iteritems():

            # Parse radar data
            data = radar.fields[field]['data']
            dist = histogram['normalized histogram']
            bins = histogram['bin centers']

    # Process radar moments
    if moments is not None:
        for label, moment in moments.iteritems():
            for field, histogram in moment.iteritems():
                continue


    return
