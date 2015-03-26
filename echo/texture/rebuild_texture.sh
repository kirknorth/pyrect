# Remove and rebuild _texture_fields.so

rm -f compute_texture.so
f2py -m compute_texture -h _compute_texture.pyf src/compute_texture.f90 --overwrite-signature
f2py --fcompiler=gfortran -c _compute_texture.pyf src/compute_texture.f90 
