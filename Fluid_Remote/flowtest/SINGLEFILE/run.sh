gfortran -c userfunc.f -o userfunc.o
#gfortran main.f liblapack.a librefblas.a userfunc.o mesh.o physol.o -o main
gfortran main.f liblapack.a librefblas.a userfunc.o -o main

