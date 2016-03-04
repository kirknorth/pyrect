# Remove and build _texture.so (Cython)
rm -vf _texture.so
cython _texture.pyx
python setup.py build_ext --inplace

# Remove and build _texture.so (Fortran)
rm -vf _texture.so
f2py -m _texture -h _texture.pyf _texture.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _texture.pyf _texture.f90 