"""
"""

from os.path import join
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):
    """
    """
    config = Configuration('texture', parent_package, top_path)

    # Add texture fields extension
    source = ['_compute_texture.pyf', 'src/compute_texture.f90']
    config.add_extension(
        'compute_texture', sources=source, f2py_options=[])

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
