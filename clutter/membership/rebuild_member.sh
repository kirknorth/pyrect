# Remove and rebuild _texture_fields.so

rm -f member.so
f2py -m member -h _member.pyf src/member.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _member.pyf src/member.f90 
