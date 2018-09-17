
CC = g++

FLAGS = -std=c++11 -fopenmp  -O3
INCLUDE = -I$(HOME)/include -I$(HOME)/include/armadillo -I$(HOME)/imsrg/work/compiled -I$(HOME)/imsrg/src
#LIBS = -lIMSRG -lopenblas -llapack -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS = $(HOME)/imsrg/src/libIMSRG.so -lopenblas  -lgsl -lz -lhdf5 -lhdf5_cpp
#LIBS = -L$(HOME)/lib64 -L$(HOME)/lib -lopenblas -lgsl -lz -lhdf5 -lhdf5_cpp
LIBS += -llapack
ifeq ($(DEBUG),on)
 FLAGS += -g
endif

all: exe
exe: Atomic writeAtomicTBME

%:	%.cc
	$(CC) $(INCLUDE) $^ -o $@ $(FLAGS) $(LIBS)

clean:
	rm -f *.o *.so boost_src/*.o
