# Remove and rebuild _texture_fields.so

rm _texture_fields.so
f2py -m texture_fields -h _texture_fields.pyf src/texture_fields.f90
f2py --fcompiler=gfortran -c _texture_fields.pyf src/texture_fields.f90 