gcc -c -g cube_test.c -fPIC 
gcc -o cube_test cube_test.o -lgsl cube/hcubature.c  -lgslcblas -lm -g -fPIC
#./cube_test 
