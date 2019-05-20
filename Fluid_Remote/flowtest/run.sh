gfortran -c userfunc.f -o userfunc.o
gfortran -c mesh.f -o mesh.o
gfortran -c physol.f -o physol.o
#gfortran main.f liblapack.a librefblas.a userfunc.o mesh.o physol.o -o main
gfortran main.f userfunc.o mesh.o physol.o -o main

