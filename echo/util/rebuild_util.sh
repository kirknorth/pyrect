# Remove and rebuild _texture_fields.so

rm -f util_brute.so
f2py -m util_brute -h _util_brute.pyf src/util_brute.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _util_brute.pyf src/util_brute.f90 
