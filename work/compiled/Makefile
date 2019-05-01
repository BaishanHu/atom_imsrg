
CC = g++

FLAGS = -std=c++11 -fopenmp  -O3
INCLUDE = -I$(HOME)/include -I$(HOME)/include/armadillo -I./compiled -I../../src
#LIBS = -lIMSRG -lopenblas -llapack -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS = $(HOME)/lib/libIMSRG.so -lopenblas  -lgsl -lz -lhdf5 -lhdf5_cpp
#LIBS = -L$(HOME)/lib64 -L$(HOME)/lib -lopenblas -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS += -lmkl_sequential -lmkl_core -lgomp
ifeq ($(DEBUG),on)
 FLAGS += -g
endif

all: exe
exe: Atomic writeAtomicTBME

%:	%.cc
	$(CC) $(INCLUDE) $^ -o $@ $(FLAGS) $(LIBS)

clean:
	rm -f *.o *.so boost_src/*.o
