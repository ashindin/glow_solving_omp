ifort -O2 -c glow_solving_omp_test.f90
ifort -O2 -c sphere_routs.f90 asa047_sphere.f90 fitsio.f fitsfort.f fitsf90.f
ifort -O2 -o glow_solving_omp_test.out glow_solving_omp_test.o sphere_routs.o asa047_sphere.o fitsio.o fitsfort.o fitsf90.o -qopenmp

