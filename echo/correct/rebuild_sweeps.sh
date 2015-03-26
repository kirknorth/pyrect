# Remove and rebuild _texture_fields.so

rm -f sweeps.so
f2py -m sweeps -h _sweeps.pyf src/sweeps.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _sweeps.pyf src/sweeps.f90 
