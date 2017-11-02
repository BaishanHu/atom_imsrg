gcc -c -g qagiu_test.c
gcc -o qagiu_test qagiu_test.o -lgsl -lgslcblas -lm -g
./qagiu_test 

#gcc -O3 -g -std=c11 -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wold-style-definition -lgsl -lgslcblas -lm qag.c -o qag
