# Remove and rebuild _texture_fields.so

rm -f _compute_texture.so
f2py -m compute_texture -h _compute_texture.pyf src/compute_texture.f90
f2py --fcompiler=gfortran -c _compute_texture.pyf src/compute_texture.f90 
