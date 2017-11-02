gcc -c -g mc_test.c
gcc -o mc_test mc_test.o -lgsl -lgslcblas -lm -g
./mc_test 
