# Remove and build _noise.so
rm -vf _noise.so
f2py -m _noise -h _noise.pyf _noise.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _noise.pyf _noise.f90 
