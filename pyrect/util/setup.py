"""
"""

import os
import sys
from numpy import get_include
from numpy.distutils.misc_util import Configuration
from numpy.distutils.core import setup


def configuration(parent_package='', top_path=None):
    """
    """

    config = Configuration('util', parent_package, top_path)

    # Add extension for _texture submodule (Cython)
    config.add_extension(
        '_texture', sources='_texture.c', include_dirs=get_include())

    # Add extension for _texture submodule (Fortran)
    config.add_extension(
        '_texture', sources=['_texture.pyf', '_texture.f90'],
        f2py_options=None, include_dirs=get_include())

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
