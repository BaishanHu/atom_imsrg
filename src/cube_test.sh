g++ -c -g cube_test.cc -fPIC 
g++ -o cube_test cube_test.o -lgsl -lgslcblas -lm -g -fPIC
./cube_test 
