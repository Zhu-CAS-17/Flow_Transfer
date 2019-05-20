gfortran -c userfunc.f -o userfunc.o
gfortran testmatrix.f userfunc.o -o testmatrix
./testmatrix
