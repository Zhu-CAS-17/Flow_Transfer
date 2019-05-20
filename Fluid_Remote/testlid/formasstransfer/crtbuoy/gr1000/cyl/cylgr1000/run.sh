gfortran -c userfunc.f -o userfunc.o
#gfortran -c phyproj.f -llibfishpack.a -o phyproj.o
#gfortran main.f liblapack.a librefblas.a userfunc.o mesh.o physol.o -o main
gfortran main.f libfish90.a userfunc.o -o cylgr1000
#gfortran main.f userfunc.o phyproj.o -o main

