gcc -c -g radIntTest.c 
gcc -o radIntTest radIntTest.o -lgsl -lgslcblas -lm
./radIntTest
