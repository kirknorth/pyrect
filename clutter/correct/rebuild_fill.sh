# Remove and rebuild _texture_fields.so

rm -f fill_holes.so
f2py -m fill_holes -h _fill_holes.pyf src/fill_holes.f90
f2py --fcompiler=gfortran -c _fill_holes.pyf src/fill_holes.f90 
